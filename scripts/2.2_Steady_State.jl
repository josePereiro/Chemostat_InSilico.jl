import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement
    const InU = InCh.Utilities

    using ProgressMeter
    using Plots
    using Plots.PlotMeasures
    using InteractiveUtils

    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Serialization
    using Statistics
    using Dates
end

## ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
fig_path(fname) = joinpath(InLP.DYN_FIGURES_DIR, fname)
minmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

## ----------------------------------------------------------------------------
# Load and clear DAT
# DAT [[:DyM, :TS, :M], Vl, D, ϵ]
DAT = UJL.DictTree()
let

    RAWDAT_FILE = InCh.DYN_DATA_BUNDLE_FILE
    RAWDAT = UJL.load_data(RAWDAT_FILE)
    # RAWDAT[:death_th] = 1e-2
    # UJL.save_data(RAWDAT_FILE, RAWDAT)

    # Clear deaths
    for Vl in RAWDAT[:Vls], D in RAWDAT[:Ds], ϵ in RAWDAT[:ϵs]
        
        M, TS = RAWDAT[[:M, :TS], Vl, D, ϵ]
        
        # Check dead
        M.X < RAWDAT[:death_th] && continue

        DAT[:TS, Vl, D, ϵ] = TS
        DAT[:M0, Vl, D, ϵ] = M

        push!(get!(DAT, [], :Ds), M.D)
        push!(get!(DAT, [], :Vls), M.Vl)
        push!(get!(DAT, [], :ϵs), M.ϵ)
    end
    
    sort!(unique!(DAT[:Ds]))
    sort!(unique!(DAT[:Vls]))
    sort!(unique!(DAT[:ϵs]))

end

## ----------------------------------------------------------------------------
# Prepare network
const STST_POL = :STST_POL   # Use polytope computed from chemostat stst assertion
const DYN_POL = :DYN_POL     # Use dynamic polytope
const HOMO = :HOMO           # Do not use extra constraints
const HETER = :HETER         # Match ME and Dy biom average
const FIXXED = :FIXXED       # Fix biom around observed

MOD_COLORS = Dict(
    HOMO => :red, 
    FIXXED => :blue,
    HETER => :orange,
)

POL_STYLE = Dict(
    STST_POL => :dot,
    DYN_POL => :dash,
)

function prepare_pol!(M, POLTsym)
    net = M.net
    if POLTsym == DYN_POL
        # Use Michaelis-Menden
        net.ub[M.vg_idx] = max(net.lb[M.vg_idx], (M.Vg * M.sg) / (M.Kg + M.sg))
        net.ub[M.vl_idx] = max(net.lb[M.vl_idx], (M.Vl * M.sl) / (M.Kl + M.sl))
    elseif  POLTsym == STST_POL
        # Use Chemostat Steady-State assumption
        # ub < max(V, c/ xi) see cossio's paper
        xi = M.X / M.D
        net.ub[M.vg_idx] = max(net.lb[M.vg_idx], min(M.Vg , M.cg / xi))
        net.ub[M.vl_idx] = max(net.lb[M.vl_idx], min(M.Vl , M.cl / xi))
    else
        error("Unknown pol type")
    end
    return net
end

function run_model!(M, MODsym; LP_cache, δ, δμ, DyBiom)
    
    # MaxEnt fun
    y = InLP.Y # atp/biomass yield
    maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)

    if MODsym == HETER
        # Gradient descent
        target = DyBiom
        x0 = 5e2
        x1 = x0 * 0.9
        maxΔ = x1 * 0.5
        th = 1e-3
        maxiters = 500
        
        # Dynamic caching
        beta0 = InU.grad_desc(;target, x0, x1, maxΔ, th, maxiters) do beta
            MEMs = InLP.get_marginals(maxentf(beta), M, [InLP.BIOMASS_IDER]; 
                δ, verbose = false, LP_cache)
            return InLP.av(MEMs[InLP.BIOMASS_IDER])
        end
    else
        beta0 = 0.0
    end

    if MODsym == FIXXED
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = DyBiom * (1.0 + δμ)
        net.lb[M.obj_idx] = DyBiom * (1.0 - δμ)
    end

    MEMs = InLP.get_marginals(maxentf(beta0), M; δ, LP_cache)
    return MEMs, beta0
end

## ----------------------------------------------------------------------------
# COMPUTE MARGINALS
let

    δ = DAT[:δ]   = 0.08 # marginal discretization factor
    δμ = DAT[:δμ] = 0.001 # FIXXED biomass variance

    N = prod(length.(DAT[[:Vls, :Ds, :ϵs]]))
    c = 0
    
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]
        c += 1

        !haskey(DAT, :M0, Vl, D, ϵ) && continue
        M0, TS = DAT[[:M0, :TS], Vl, D, ϵ]
        LP_cache = InLP.vgvatp_cache(M0)
        @info "Doing $c/$N ... " Vl D ϵ M0.X

        ## ----------------------------------------------------------------------------
        # Dynamic marginal
        f(vatp, vg) = M0.Xb[vatp][vg] / M0.X
        DyMs = DAT[:DyMs, Vl, D, ϵ] = InLP.get_marginals(f, M0; δ, LP_cache)
        DyBiom = InLP.av(DyMs[InLP.BIOMASS_IDER]) # biomass dynamic mean

        ## ----------------------------------------------------------------------------
        # MaxEnt marginals
        for POLTsym in [STST_POL, DYN_POL]
            for MODsym in [HOMO, FIXXED, HETER]
                
                # Setup network
                M = deepcopy(M0)
                prepare_pol!(M, POLTsym)
                
                MEMs, beta0 = run_model!(M, MODsym; LP_cache, δ, δμ, DyBiom)
                
                DAT[MODsym, POLTsym, :M, Vl, D, ϵ] = M
                DAT[MODsym, POLTsym, :MEMs, Vl, D, ϵ] = MEMs
                DAT[MODsym, POLTsym, :beta0, Vl, D, ϵ] = beta0

                # Ranges
                vatp_range, vg_ranges = InLP.vatpvg_ranges(M)
                DAT[MODsym, POLTsym, :POL, Vl, D, ϵ] = (;vatp_range, vg_ranges)
            end
        end
        println()
    end
end


## ----------------------------------------------------------------------------
# PLOTS
PS = UJL.DictTree()

# ## ----------------------------------------------------------------------------
# # Bound Correlation
# let
#     f(x) = x
#     p = plot(;title = "Bounds Correlation", 
#         xlabel = "dym bound", ylabel = "stst bound")
#     l, u = Inf, -Inf
#     for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

#         !haskey(DAT, :M0, Vl, D, ϵ) && continue
#         M0 = DAT[:M, Vl, D, ϵ]
#         xi = M0.X / D

#         vgubs = []
#         vlubs = []
#         for (MEMs_sym, prep_fun!) in [
#                                     (:MEMs_dym_bound, prepare_net_dym!),
#                                     (:MEMs_stst_bound, prepare_net_stst!),
#                                 ]

#             M = deepcopy(M0)
#             InLP.fill_board!(M.Xb, 1.0)
#             net = prep_fun!(M, xi)
#             push!(vgubs, net.ub[M.vg_idx])
#             push!(vlubs, net.ub[M.vl_idx])

#         end
#         xs = first.([vgubs, vlubs])
#         ys = last.([vgubs, vlubs])
#         l = minimum([l; xs; ys])            
#         u = maximum([u; xs; ys])            
#         scatter!(p, f.(xs), f.(ys); alpha = 0.5, color = :black, label = "")

#     end
#     plot!(p, f.([l, u]), f.([l, u]); label = "", ls = :dash, alpha = 0.8)

#     # saving
#     pname = UJL.mysavename("bounds_corr", "png")
#     P[pname] = deepcopy(p)
#     fname = fig_path(string(fileid, "_", pname))    
#     savefig(p, fname)
#     @info "Plotting" fname

# end

## ----------------------------------------------------------------------------
function plot_pol(POLTsym, MODsym, Vl, D, ϵ; gparams = Dict(), sparams = Dict())
    
    p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg", gparams...)
    # Polytopes
    vatp_range, vg_ranges = DAT[MODsym, POLTsym, :POL, Vl, D, ϵ] 
    vatps, vgLs, vgUs = [], [], []
    
    for (vatpi, vatp) in enumerate(vatp_range)
        vg_range = vg_ranges[vatpi]
        isempty(vg_range) && continue
        vgL, vgU = minmax(vg_range)
        push!(vatps, vatp)
        push!(vgLs, vgL)
        push!(vgUs, vgU)
    end

    color = MOD_COLORS[MODsym]
    plot!(p, [vatps], [vgLs]; label = "", color, lw = 3, sparams...)
    plot!(p, [vatps], [vgUs]; label = "", color, lw = 3, sparams...)
    PS[MODsym, POLTsym, :POL, Vl, D, ϵ] = deepcopy(p)
end


## ----------------------------------------------------------------------------
function plot_marginals(MODsym, POLTsym, Vl, D, ϵ)

    M = DAT[MODsym, POLTsym, :M, Vl, D, ϵ]
    DyMs = DAT[:DyMs, Vl, D, ϵ]

    sparams =(;alpha = 0.8, lw = 3, ylim = [0.0, Inf])
    gparams = (xaxis = nothing, yaxis = nothing, grid = false, 
            titlefont = 10, xaxisfont = 10)

    ps = UJL.DictTree()
    for rxn in M.net.rxns
        p = ps[rxn] = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
        plot!(p, DyMs[rxn]; label = "", color = :black, sparams...)
    end

    ls = :solid #POL_STYLE[POLTsym]
    color = MOD_COLORS[MODsym]
    MEMs = DAT[MODsym, POLTsym, :MEMs, Vl, D, ϵ]

    # Marginals
    for rxn in M.net.rxns
        try; plot!(ps[rxn], MEMs[rxn]; label = "", color, sparams...)
            catch err; @error string("Doing $rxn\n", UJL.err_str(err))
        end
    end

    # Polytopes
    ps[:pol] = plot_pol(POLTsym, MODsym, Vl, D, ϵ; gparams, sparams)

    vps = [ps[M.net.rxns]; ps[:pol]]
    p = plot(vps...; layout = length(vps))

    # saving
    pname = UJL.mysavename("Dyn_marginals_$(POLTsym)_$(MODsym)", "png"; Vl, D, ϵ)
    fname = fig_path(string(fileid, "_", pname))    
    savefig(p, fname)
    @info "Plotting" fname

end
let
    # Vl = DAT[:Vls] |> first
    # D = DAT[:Ds][2]
    # ϵ = DAT[:ϵs] |> first
    # MODsym = HETER
    # POLTsym = STST_POL
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]
        for POLTsym in [STST_POL, DYN_POL], MODsym in [HOMO, FIXXED, HETER]
            plot_marginals(MODsym, POLTsym, Vl, D, ϵ)
        end
    end
end

## ----------------------------------------------------------------------------
# plot marginals
let
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]
        
        
        !haskey(DAT, :M0, Vl, D, ϵ) && continue
        M = DAT[:M, Vl, D, ϵ]
        DyMs = DAT[:DyMs, Vl, D, ϵ]
        
        sparams =(;alpha = 0.8, lw = 3, ylim = [0.0, Inf])
        gparams = (xaxis = nothing, yaxis = nothing, grid = false, 
        titlefont = 10, xaxisfont = 10)
        
        ps = UJL.DictTree()
        ps[:pol] = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg", gparams...)
        for rxn in M.net.rxns
            p = ps[rxn] = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            plot!(p, DyMs[rxn]; label = "", color = :black, sparams...)
        end


        for POLTsym in [STST_POL]
            ls = :solid # POL_STYLE[POLTsym]
            for MODsym in [HOMO, FIXXED, HETER]

                MEMs = DAT[MODsym, POLTsym, Vl, D, ϵ]
                color = MOD_COLORS[MODsym]

                # Marginals
                for rxn in M.net.rxns
                    try; plot!(ps[rxn], MEMs[rxn]; label = "", color, ls, sparams...)
                        catch err; @error string("Doing $rxn\n", UJL.err_str(err))
                    end
                end

                # Polytopes
                vatp_range, vg_ranges = DAT[MODsym, POLTsym, :ranges, Vl, D, ϵ] 
                vatps, vgLs, vgUs = [], [], []
                for (vatpi, vatp) in enumerate(vatp_range)
                    vg_range = vg_ranges[vatpi]
                    isempty(vg_range) && continue
                    vgL, vgU = minmax(vg_range)
                    push!(vatps, vatp)
                    push!(vgLs, vgL)
                    push!(vgUs, vgU)
                end
                plot!(ps[:pol], [vatps], [vgLs]; label = "", color, sparams...)
                plot!(ps[:pol], [vatps], [vgUs]; label = "", color, sparams...)

            end
        end # for POLTsym


        ps = [[ps[rxn] for rxn in M.net.rxns]; ps[:pol]]
        p = plot(ps...; layout = length(ps))

        # saving
        pname = UJL.mysavename("Dyn_marginals", "png"; Vl, D, ϵ)
        fname = fig_path(string(fileid, "_", pname))    
        savefig(p, fname)
        @info "Plotting" fname

    end # for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]
end

## ----------------------------------------------------------------------------
# Steady State EP Dynamic correlation
let
    f(x) = x

    for MEMs_sym in [:MEMs_dym_bound, :MEMs_stst_bound]

        p = plot(;title = "Total Correlation", 
            xlabel = "flx (Dynamic)", ylabel = "flx (MaxEnt)")

        PD = UJL.DictTree()
        @showprogress for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

            !haskey(DAT, :M0, Vl, D, ϵ) && continue
            M, TS = DAT[[:M, :TS], Vl, D, ϵ]
            DyMs = DAT[:DyMs, Vl, D, ϵ]
            MEMs = DAT[MEMs_sym, Vl, D, ϵ]
            beta0 = DAT[MEMs_sym, :beta0, Vl, D, ϵ]
            isnan(beta0) && continue # leave out deaths

            for rxn in M.net.rxns
                rxn in ["atpm", "vatp"] && continue

                foreach([:MEAv, :MEStd, :DyAv, :DyStd]) do id 
                    get!(PD, [], rxn, id)
                end

                push!(PD[rxn, :MEAv], InLP.av(MEMs[rxn]))
                push!(PD[rxn, :MEStd], sqrt(InLP.va(MEMs[rxn])))
                push!(PD[rxn, :DyAv], InLP.av(DyMs[rxn]))
                push!(PD[rxn, :DyStd], sqrt(InLP.va(DyMs[rxn])))

            end
        end

        # scatter!(p, f.([DyAv]), f.([MEAv]); xerr = f.([DyStd]), yerr = f.([MEStd]),
        #             alpha = 0.5, color = :white, label = "", ms = 8)

        for rxn in keys(PD)
            norm = maximum(abs.(PD[rxn, :DyAv]))
            scatter!(p, PD[rxn, :DyAv] ./ norm, PD[rxn, :MEAv] ./ norm; 
                    xerr = PD[rxn, :DyStd] ./ norm, yerr = PD[rxn, :MEStd] ./ norm,
                    alpha = 0.5, label = "", ms = 8)
        end

        plot!(p, [-1.0, 1.0], [-1.0, 1.0]; label = "", ls = :dash, alpha = 0.8)

        # saving
        pname = UJL.mysavename("Dyn_$(MEMs_sym)_flx_corr", "png")
        fname = fig_path(string(fileid, "_", pname))    
        savefig(p, fname)
        @info "Plotting" fname
    end
end