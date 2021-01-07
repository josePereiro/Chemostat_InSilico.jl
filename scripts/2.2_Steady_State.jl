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
minmax(a::Vector) = (minimum(a), maximum(a))

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

# ----------------------------------------------------------------------------
# Globals
δ = DAT[:δ] = 0.08 # marginal discretization factor

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
    net.ub[M.vg_idx] = max(net.lb[M.vg_idx], max(M.Vg , M.cg / xi))
    net.ub[M.vl_idx] = max(net.lb[M.vl_idx], max(M.Vl , M.cl / xi))
    net
end

## ----------------------------------------------------------------------------
# TO DELETE, fix UtilsJL dep
# Base.get(pd::UJL.DictTree, defl, k, ks...) = haskey(pd, k, ks...) ? pd[k, ks...] : defl

## ----------------------------------------------------------------------------
# COMPUTE MARGINALS
let
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

        M0, TS = DAT[[:M, :TS], Vl, D, ϵ]
        LP_cache = InLP.vgvatp_cache(M0)
        @info "Doing" Vl D ϵ M0.X

        ## ----------------------------------------------------------------------------
        # Dynamic marginal
        f(vatp, vg) = M0.Xb[vatp][vg] / M0.X
        DyMs = DAT[:DyMs, Vl, D, ϵ] = InLP.get_marginals(f, M0; δ, LP_cache)

        ## ----------------------------------------------------------------------------
        # MaxEnt marginals
        for (prep_mod, prep_fun!) in [
                                    (:dym, prepare_net_dym!),
                                    (:stst, prepare_net_dym!),
                                ]

            MEM_sym = Symbol(string("MEM", prep_mod))
            
            M = deepcopy(M0)
            xi = M.X / D

            # Setup network
            net = prep_fun!(M, xi)

            # MaxEnt fun
            y = InLP.Y # atp/biomass yield
            maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)
            
            # # Gradient descent
            biom_ider = InLP.BIOMASS_IDER
            target = InLP.av(DyMs[biom_ider]) # biomass dynamic mean
            x0 = get(DAT, 1e3, MEM_sym, :beta0, Vl, D, ϵ) 
            x1 = x0 * 1.1
            maxΔ = x0 * 0.5
            th = 1e-3
            maxiters = 50
            # Dynamic caching
            beta0 = InU.grad_desc(;target, x0, x1, maxΔ, th, maxiters) do beta
                MEMs = InLP.get_marginals(maxentf(beta), M, [biom_ider]; 
                    δ, verbose = false, LP_cache)
                f = InLP.av(MEMs[biom_ider])
            end
            DAT[MEM_sym, :beta0, Vl, D, ϵ] = beta0
            DAT[MEM_sym, Vl, D, ϵ] = InLP.get_marginals(maxentf(beta0), M; δ, LP_cache)
        end
            
    end
end


#  ----------------------------------------------------------------------------
# PLOTS

## ----------------------------------------------------------------------------
# state_vs_D
let
    f(x) = x

    for Vl in DAT[:Vls], ϵ in DAT[:ϵs]


        p = plot(;title = "", 
            xlabel = "flx (Dynamic)", ylabel = "flx (MaxEnt)")

        DT = UJL.DictTree()
        for D in DAT[:Ds]

            M = DAT[:M, Vl, D, ϵ]
            DyMs = DAT[:DyM, Vl, D, ϵ]
            MEMs = DAT[:MEMs, Vl, D, ϵ]
            beta0 = DAT[:beta0, Vl, D, ϵ]
            isnan(beta0) && continue # leave out deaths

            for rxn in M.net.rxns
                #
            end

        end

        # saving
        pname = UJL.mysavename("Dy_ME_total_correlation", "png")
        fname = fig_path(string(fileid, "_", pname))    
        savefig(p, fname)
        @info "Plotting" fname
    end
end

## ----------------------------------------------------------------------------
# plot marginals
let
    for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

        M, TS = DAT[[:M, :TS], Vl, D, ϵ]
        DyMs = DAT[:DyM, Vl, D, ϵ]
        MEMs = DAT[:MEMs, Vl, D, ϵ]
        beta0 = DAT[:beta0, Vl, D, ϵ]
        isnan(beta0) && continue # leave out deaths

        ps = []
        gparams = (;xlabel = "flx", ylabel = "prob", 
            xaxis = nothing, yaxis = nothing, grid = false, 
            titlefont = 10, xaxisfont = 10)
        sparams =(;alpha = 0.8, lw = 3, ylim = [0.0, Inf])
        for rxn in M.net.rxns
            DyM = DyMs[rxn]
            MEM = MEMs[rxn]
            
            try
                p = plot(;title = rxn, gparams...)
                plot!(p, DyM; label = "", color = :red, sparams...)
                plot!(p, MEM; label = "", color = :blue, sparams...)
                push!(ps, p)
            catch err 
                @error string("Doing $rxn\n", UJL.err_str(err))
            end
        end
        p = plot(ps...; layout = length(ps))

        # saving
        pname = UJL.mysavename("marginals", "png"; Vl, D, ϵ, beta0)
        fname = fig_path(string(fileid, "_", pname))    
        savefig(p, fname)
        @info "Plotting" fname

    end
end

## ----------------------------------------------------------------------------
# Steady State EP Dynamic correlation
let
    f(x) = x
    p = plot(;title = "Total Correlation", 
        xlabel = "flx (Dynamic)", ylabel = "flx (MaxEnt)")
    l, u = Inf, -Inf
    @showprogress for Vl in DAT[:Vls], D in DAT[:Ds], ϵ in DAT[:ϵs]

        M, TS = DAT[[:M, :TS], Vl, D, ϵ]
        DyMs = DAT[:DyM, Vl, D, ϵ]
        MEMs = DAT[:MEMs, Vl, D, ϵ]
        beta0 = DAT[:beta0, Vl, D, ϵ]
        isnan(beta0) && continue # leave out deaths

        for rxn in M.net.rxns
            rxn == "atpm" && continue
            DyAv = InLP.av(DyMs[rxn])
            DyStd = sqrt(InLP.va(DyMs[rxn]))
            MEAv = InLP.av(MEMs[rxn])
            MEStd = sqrt(InLP.va(MEMs[rxn]))

            scatter!(p, f.([DyAv]), f.([MEAv]); xerr = f.([DyStd]), yerr = f.([MEStd]),
                alpha = 0.5, color = :black, label = "")

            l = minimum([l, DyAv, MEAv])            
            u = maximum([u, DyAv, MEAv])            
            isnan(l) && @show rxn, l, DyAv, MEAv
            isnan(u) && @show rxn, u, DyAv, MEAv
        end
    end

    plot!(p, f.([l, u]), f.([l, u]); label = "", ls = :dash, alpha = 0.8)

    # saving
    pname = UJL.mysavename("Dy_ME_total_correlation", "png")
    fname = fig_path(string(fileid, "_", pname))    
    savefig(p, fname)
    @info "Plotting" fname
end