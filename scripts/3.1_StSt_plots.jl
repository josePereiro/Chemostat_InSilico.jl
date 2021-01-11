import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement
    const InU = InCh.Utilities

    using ProgressMeter
    using Plots

    # Test
    # import GR
    # GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    using Base.Threads
    using Dates
    using Statistics
    using InteractiveUtils
    using Random

end

## ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
minmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

## ----------------------------------------------------------------------------
# Prepare network
const STST_POL = :STST_POL   # Use polytope computed from chemostat stst assertion
const DYN_POL = :DYN_POL     # Use dynamic polytope
const HOMO = :HOMO           # Do not use extra constraints
const HETER = :HETER         # Match ME and Dy biom average
const FIXXED = :FIXXED       # Fix biom around observed

MOD_COLORS = Dict(
    HOMO => :red, 
    FIXXED => :orange,
    HETER => :blue,
)

POL_STYLE = Dict(
    STST_POL => :dot,
    DYN_POL => :dash,
)

## ----------------------------------------------------------------------------
# Marginals
# MDAT[MODsym, POLTsym, :M, Vl, D, ϵ]
# MDAT[MODsym, POLTsym, :MEMs, Vl, D, ϵ]
# MDAT[MODsym, POLTsym, :beta0, Vl, D, ϵ]
# MDAT[:STATUS, Vl, D, ϵ]
MDAT = UJL.load_data(InCh.MARGINALS_DATA_FILE);
POLTsym = STST_POL
EXP_PARAMS = Iterators.product(MDAT[[:Vls, :Ds, :ϵs, :τs]]...);

## ----------------------------------------------------------------------------
# PLOTS
PS = UJL.DictTree()
function mysavefig(p, pname; params...)
    pname = UJL.mysavename(pname, "png"; params...)
    fname = joinpath(InLP.DYN_FIGURES_DIR, string(fileid, "_", pname))
    PS[pname] = deepcopy(p)
    savefig(p, fname)
    @info "Plotting" fname
end

## ----------------------------------------------------------------------------
# Plot functions
function plot_pol!(p, POLTsym, MODsym, Vl, D, ϵ; sparams...)
    
    vatp_range, vg_ranges = MDAT[MODsym, POLTsym, :POL, Vl, D, ϵ] 
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
    ls = POL_STYLE[POLTsym]
    plot!(p, [vatps], [vgLs]; label = "", ls, alpha = 0.4, color, lw = 3, sparams...)
    plot!(p, [vatps], [vgUs]; label = "", ls, alpha = 0.4, color, lw = 3, sparams...)
    PS[MODsym, POLTsym, :POL, Vl, D, ϵ] = deepcopy(p)
end

function plot_marginals!(p, MODsym, POLTsym, rxn, Vl, D, ϵ; sparams...)

    ls = POL_STYLE[POLTsym]
    color = MOD_COLORS[MODsym]
    MEMs = MDAT[MODsym, POLTsym, :MEMs, Vl, D, ϵ]
    DyMs = MDAT[:DyMs, Vl, D, ϵ]
    
    # Marginals
    plot!(p, DyMs[rxn]; label = "", sparams..., color = :black)
    plot!(p, MEMs[rxn]; label = "", color, sparams...)
end

## ----------------------------------------------------------------------------
# Polytopes
let
    sparams =(;alpha = 0.5, lw = 5, ylim = [0.0, Inf])

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MDAT[:STATUS, Vl, D, ϵ] == :death && continue

        for MODsym = [HOMO, FIXXED, HETER]
            p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg")        
            for POLTsym in [STST_POL, DYN_POL]
                plot_pol!(p, POLTsym, MODsym, Vl, D, ϵ; sparams...)
            end
            
            mysavefig(p, "Polytopes_$(POLTsym)_$(MODsym)"; Vl, D, ϵ)
        end
    end
end

## ----------------------------------------------------------------------------
# Marginals
let
    D = MDAT[:Ds][1]
    Vl = MDAT[:Vls][1]
    ϵ = MDAT[:ϵs][1]
    τ = MDAT[:τs][1]

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MDAT[:STATUS, Vl, D, ϵ] == :death && continue
        ps = []
        sparams =(;alpha = 0.7, lw = 5, ylim = [0.0, Inf])
        gparams = (xaxis = nothing, yaxis = nothing, grid = false, 
                titlefont = 10, xaxisfont = 10)
        for rxn in InLP.RXNS
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym = [HOMO, FIXXED, HETER]
                plot_marginals!(p, MODsym, POLTsym, rxn, Vl, D, ϵ; sparams...)
            end
            push!(ps, p)
        end

        p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg", gparams...)        
        plot_pol!(p, POLTsym, HETER, Vl, D, ϵ; sparams...)
        plot_pol!(p, POLTsym, HOMO, Vl, D, ϵ; sparams...)
        plot_pol!(p, POLTsym, FIXXED, Vl, D, ϵ; sparams...)
        push!(ps, p)
        
        p = plot(ps...; layout = length(ps))
        mysavefig(p, "Marginals_$(POLTsym)_"; Vl, D, ϵ)
    end

end


## ----------------------------------------------------------------------------
let
    D = MDAT[:Ds][1]
    Vl = MDAT[:Vls][1]
    ϵ = MDAT[:ϵs][1]
    τ = MDAT[:τs][1]
end
## ----------------------------------------------------------------------------
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

# ## ----------------------------------------------------------------------------
# # Steady State EP Dynamic correlation
# let
#     f(x) = x

#     for MEMs_sym in [:MEMs_dym_bound, :MEMs_stst_bound]

#         p = plot(;title = "Total Correlation", 
#             xlabel = "flx (Dynamic)", ylabel = "flx (MaxEnt)")

#         PD = UJL.DictTree()
#         @showprogress for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

#             !haskey(DAT, :M0, Vl, D, ϵ) && continue
#             M, TS = DAT[[:M, :TS], Vl, D, ϵ]
#             DyMs = DAT[:DyMs, Vl, D, ϵ]
#             MEMs = DAT[MEMs_sym, Vl, D, ϵ]
#             beta0 = DAT[MEMs_sym, :beta0, Vl, D, ϵ]
#             isnan(beta0) && continue # leave out deaths

#             for rxn in M.net.rxns
#                 rxn in ["atpm", "vatp"] && continue

#                 foreach([:MEAv, :MEStd, :DyAv, :DyStd]) do id 
#                     get!(PD, [], rxn, id)
#                 end

#                 push!(PD[rxn, :MEAv], InLP.av(MEMs[rxn]))
#                 push!(PD[rxn, :MEStd], sqrt(InLP.va(MEMs[rxn])))
#                 push!(PD[rxn, :DyAv], InLP.av(DyMs[rxn]))
#                 push!(PD[rxn, :DyStd], sqrt(InLP.va(DyMs[rxn])))

#             end
#         end

#         # scatter!(p, f.([DyAv]), f.([MEAv]); xerr = f.([DyStd]), yerr = f.([MEStd]),
#         #             alpha = 0.5, color = :white, label = "", ms = 8)

#         for rxn in keys(PD)
#             norm = maximum(abs.(PD[rxn, :DyAv]))
#             scatter!(p, PD[rxn, :DyAv] ./ norm, PD[rxn, :MEAv] ./ norm; 
#                     xerr = PD[rxn, :DyStd] ./ norm, yerr = PD[rxn, :MEStd] ./ norm,
#                     alpha = 0.5, label = "", ms = 8)
#         end

#         plot!(p, [-1.0, 1.0], [-1.0, 1.0]; label = "", ls = :dash, alpha = 0.8)

#         # saving
#         pname = UJL.mysavename("Dyn_$(MEMs_sym)_flx_corr", "png")
#         fname = fig_path(string(fileid, "_", pname))    
#         savefig(p, fname)
#         @info "Plotting" fname
#     end
# end