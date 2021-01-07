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
        DAT[:M, Vl, D, ϵ] = M

        push!(get!(DAT, [], :Ds), M.D)
        push!(get!(DAT, [], :Vls), M.Vl)
        push!(get!(DAT, [], :ϵs), M.ϵ)
    end
    
    sort!(unique!(DAT[:Ds]))
    sort!(unique!(DAT[:Vls]))
    sort!(unique!(DAT[:ϵs]))
end

## ----------------------------------------------------------------------------
# Modeling mode
# Use Michaelis-Menden
function prepare_net_dym!(M, bla...)
    net = M.net
    net.ub[M.vg_idx] = max(net.lb[M.vg_idx], (M.Vg * M.sg) / (M.Kg + M.sg))
    net.ub[M.vl_idx] = max(net.lb[M.vl_idx], (M.Vl * M.sl) / (M.Kl + M.sl))
    net
end

# Use Chemostat Steady-State assumption
function prepare_net_stst!(M, xi = M.X / M.D)
    net = M.net
    # ub < max(V, c/ xi) see cossio's paper
    net.ub[M.vg_idx] = max(net.lb[M.vg_idx], min(M.Vg , M.cg / xi))
    net.ub[M.vl_idx] = max(net.lb[M.vl_idx], min(M.Vl , M.cl / xi))
    net
end

## ----------------------------------------------------------------------------
# COMPUTE MARGINALS
let

    δ = DAT[:δ] = 0.08 # marginal discretization factor

    N = sum(length.(DAT[[:Vls, :Ds, :ϵs]]))
    c = 0
    
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]
        c += 1

        !haskey(DAT, :M, Vl, D, ϵ) && continue
        M0, TS = DAT[[:M, :TS], Vl, D, ϵ]
        LP_cache = InLP.vgvatp_cache(M0)
        @info "Doing $c/$N ... " Vl D ϵ M0.X

        ## ----------------------------------------------------------------------------
        # Dynamic marginal
        f(vatp, vg) = M0.Xb[vatp][vg] / M0.X
        DyMs = DAT[:DyMs, Vl, D, ϵ] = InLP.get_marginals(f, M0; δ, LP_cache)

        ## ----------------------------------------------------------------------------
        # MaxEnt marginals
        for (prep_mod, prep_fun!) in [
                                    (:dym, prepare_net_dym!),
                                    (:stst, prepare_net_stst!),
                                ]

            MEMsym = Symbol(string("MEM", prep_mod))
            
            # Setup network
            M = deepcopy(M0)
            xi = M.X / D
            prep_fun!(M, xi)

            # MaxEnt fun
            y = InLP.Y # atp/biomass yield
            maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)
            
            # # Gradient descent
            biom_ider = InLP.BIOMASS_IDER
            target = InLP.av(DyMs[biom_ider]) # biomass dynamic mean
            x1 = get(DAT, 5e2, MEMsym, :beta0, Vl, D, ϵ) 
            x0 = x1 * 0.9
            maxΔ = x1 * 0.5
            th = 1e-3
            maxiters = 500
            # Dynamic caching
            beta0 = InU.grad_desc(;target, x0, x1, maxΔ, th, maxiters) do beta
                MEMs = InLP.get_marginals(maxentf(beta), M, [biom_ider]; 
                    δ, verbose = false, LP_cache)
                f = InLP.av(MEMs[biom_ider])
            end
            DAT[MEMsym, :beta0, Vl, D, ϵ] = beta0
            DAT[MEMsym, Vl, D, ϵ] = InLP.get_marginals(maxentf(beta0), M; δ, LP_cache)

            # Ranges
            vatp_range, vg_ranges = InLP.vatpvg_ranges(M)
            DAT[MEMsym, :ranges, Vl, D, ϵ] = (;vatp_range, vg_ranges)
        end
        println()
    end
end


## ----------------------------------------------------------------------------
# PLOTS
P = Dict()

## ----------------------------------------------------------------------------
# Bound Correlation
let
    f(x) = x
    p = plot(;title = "Bounds Correlation", xlabel = "dym bound", ylabel = "stst bound")
    l, u = Inf, -Inf
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

        !haskey(DAT, :M, Vl, D, ϵ) && continue
        M0 = DAT[:M, Vl, D, ϵ]
        xi = M0.X / D

        vgubs = []
        vlubs = []
        for (MEMsym, prep_fun!) in [
                                    (:MEMdym, prepare_net_dym!),
                                    (:MEMstst, prepare_net_stst!),
                                ]

            M = deepcopy(M0)
            InLP.fill_board!(M.Xb, 1.0)
            net = prep_fun!(M, xi)
            push!(vgubs, net.ub[M.vg_idx])
            push!(vlubs, net.ub[M.vl_idx])

        end
        xs = first.([vgubs, vlubs])
        ys = last.([vgubs, vlubs])
        l = minimum([l; xs; ys])            
        u = maximum([u; xs; ys])            
        scatter!(p, f.(xs), f.(ys); alpha = 0.5, color = :black, label = "")

    end
    plot!(p, f.([l, u]), f.([l, u]); label = "", ls = :dash, alpha = 0.8)

    # saving
    pname = UJL.mysavename("bounds_corr", "png")
    P[pname] = deepcopy(p)
    fname = fig_path(string(fileid, "_", pname))    
    savefig(p, fname)
    @info "Plotting" fname

end

## ----------------------------------------------------------------------------
# state_vs_D
let
    f(x) = x
    
    for Vl in DAT[:Vls], ϵ in DAT[:ϵs]

        p = plot(;title = "", 
            xlabel = "D", ylabel = "X")

        Ds = filter((D) -> haskey(DAT, :M, Vl, D, ϵ), DAT[:Ds])
        Ms = DAT[:M, Vl, Ds, ϵ]

        # Dynamic
        plot!(p, getfield.(Ms, :D), getfield.(Ms, :X), label = "")
            
        # MaxEnt
        for MEMsym in [:MEMdym, :MEMstst]

            xis = getfield.(Ms, :X) ./ Ds
            MEMss = DAT[MEMsym, Vl, Ds, ϵ]
            
            PD = UJL.DictTree()
            for D in DAT[:Ds]

                !haskey(DAT, :M, Vl, D, ϵ) && continue
                M = DAT[:M, Vl, D, ϵ]
                DyMs = DAT[:DyMs, Vl, D, ϵ]
                MEMs = DAT[MEMsym, Vl, D, ϵ]

                for rxn in M.net.rxns
                    
                end

            end

        end

        # saving
        pname = UJL.mysavename("Steady_State_X_vs_D", "png"; Vl, ϵ)
        fname = fig_path(string(fileid, "_", pname))    
        savefig(p, fname)
        @info "Plotting" fname
    end
end


## ----------------------------------------------------------------------------
# plot marginals
let
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]
        
        ps = Dict()

        !haskey(DAT, :M, Vl, D, ϵ) && continue
        M = DAT[:M, Vl, D, ϵ]
        DyMs = DAT[:DyMs, Vl, D, ϵ]

        sparams =(;alpha = 0.8, lw = 3, ylim = [0.0, Inf])
        gparams = (xaxis = nothing, yaxis = nothing, grid = false, 
                titlefont = 10, xaxisfont = 10)

        for rxn in M.net.rxns
            p = ps[rxn] = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            plot!(p, DyMs[rxn]; label = "", color = :black, sparams...)
        end

        ps[:pol] = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg", gparams...)
        for (MEMsym, color) in [(:MEMdym, :red), (:MEMstst, :blue)]

            MEMs = DAT[MEMsym, Vl, D, ϵ]
            # beta0 = DAT[MEMsym, :beta0, Vl, D, ϵ]

            # Marginals
            for rxn in M.net.rxns
                try
                    plot!(ps[rxn], MEMs[rxn]; label = "", color = color, sparams...)
                catch err 
                    @error string("Doing $rxn\n", UJL.err_str(err))
                end
            end

            # Polytopes
            vatp_range, vg_ranges = DAT[MEMsym, :ranges, Vl, D, ϵ] 
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

        end # for (MEMsym, color)


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

    for MEMsym in [:MEMdym, :MEMstst]

        p = plot(;title = "Total Correlation", 
            xlabel = "flx (Dynamic)", ylabel = "flx (MaxEnt)")
        
        # l, u = Inf, -Inf

        PD = UJL.DictTree()
        @showprogress for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

            !haskey(DAT, :M, Vl, D, ϵ) && continue
            M, TS = DAT[[:M, :TS], Vl, D, ϵ]
            DyMs = DAT[:DyMs, Vl, D, ϵ]
            MEMs = DAT[MEMsym, Vl, D, ϵ]
            beta0 = DAT[MEMsym, :beta0, Vl, D, ϵ]
            isnan(beta0) && continue # leave out deaths

            for rxn in M.net.rxns
                rxn in ["atpm", "vatp"] && continue
                DyAv = InLP.av(DyMs[rxn])
                DyStd = sqrt(InLP.va(DyMs[rxn]))
                MEAv = InLP.av(MEMs[rxn])
                MEStd = sqrt(InLP.va(MEMs[rxn]))

                foreach([:MEAv, :MEStd, :DyAv, :DyStd]) do id 
                    get!(PD, [], rxn, id)
                end

                push!(PD[rxn, :MEAv], MEAv)
                push!(PD[rxn, :MEStd], MEStd)
                push!(PD[rxn, :DyAv], DyAv)
                push!(PD[rxn, :DyStd], DyStd)

                # l = minimum([l, DyAv, MEAv])            
                # u = maximum([u, DyAv, MEAv])            
            end
        end

        # scatter!(p, f.([DyAv]), f.([MEAv]); xerr = f.([DyStd]), yerr = f.([MEStd]),
        #             alpha = 0.5, color = :white, label = "", ms = 8)

        for rxn in keys(PD)
            norm = maximum(abs.(PD[rxn, :DyAv]))
            scatter!(p, PD[rxn, :DyAv] ./ norm, PD[rxn, :MEAv] ./ norm; 
                    # xerr = PD[rxn, :DyStd] ./ norm, yerr = PD[rxn, :MEStd] ./ norm,
                    alpha = 0.5, label = "", ms = 8)
        end

        # plot!(p, f.([l, u]), f.([l, u]); label = "", ls = :dash, alpha = 0.8)
        plot!(p, [-1.0, 1.0], [-1.0, 1.0]; label = "", ls = :dash, alpha = 0.8)

        # saving
        pname = UJL.mysavename("Dyn_$(MEMsym)_flx_corr", "png")
        fname = fig_path(string(fileid, "_", pname))    
        savefig(p, fname)
        @info "Plotting" fname
    end
end