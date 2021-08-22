@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        SimD3, run_simD3!, check_stst, hist, 
        n_ave_conv!, tail_ave,
        Container, vec!, Vi, 
        normalizeP!, 
        save_simdat, load_simdat, simdat_file,
        set_status, get_status,
        load_batch, save_batch,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS, collect_ts
        
    import UtilsJL
    const PltU = UtilsJL.PlotsUtils
    const Ass = UtilsJL.ProjAssistant
    const GU = UtilsJL.GeneralUtils

    using Plots
    import GR
    !isinteractive() && GR.inline("png")

    using Base.Threads
    using ProgressMeter
    using BenchmarkTools
    using ExtractMacro

    import DataFileNames
    const DFN = DataFileNames

end

## ------------------------------------------------------
# globlas
# batch = (; push_frec, Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts)
Ds, ϵs, cgs, simid = Dyn.lglob(:Ds, :ϵs, :cgs, :SimD3Id)

## ------------------------------------------------------
let
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    iter = collect(Iterators.product(Ds, ϵs, cgs))
    iter = filter(iter) do tp
        D, ϵ, cg = tp
        isinf(cg) && return false
        status = Dyn.get_status(simid, (;D, ϵ, cg))
        (status != STST_SIM_STATUS) && return false
        return true
    end

    for (simi, (D, ϵ, cg)) in collect(enumerate(iter))
        
        # # Test
        # D=3.86e-01 
        # cg=1.50e+01 
        # ϵ=1.00e-02
        
        simparams = (;D, ϵ, cg)
        !isfile(Dyn.maxent_file(simid, MEmode, simparams)) && continue

        @info("At", simi)

        dat = Dyn.load_maxent(simid, MEmode, simparams)
        PME, z_beta, ug_beta, z_avPME, ug_avPME = dat

        S = Dyn.load_simdat(simid, simparams)
        PX = S.P
        V = S.V
        Vz = Vi(V, :z)
        Vug = Vi(V, :ug)
        Vuo = Vi(V, :uo)

        # Plots
        ps = Plots.Plot[]
        for (title, P) in [("Dyn", PX), ("MaxEnt", PME)]
            Pz = hist(Vz, P)
            p = bar(Pz; label = "", title,
                xlabel = "z", ylabel = "proj"
            )
            z_av = sum(P .* Vz)
            vline!(p, [z_av]; label = "", color = :red)
            push!(ps, p)
            
            Pug = hist(Vug, P)
            p = bar(Pug; label = "", title,
                xlabel = "ug", ylabel = "proj"
            )
            ug_av = sum(P .* Vug)
            vline!(p, [ug_av]; label = "", color = :red)
            push!(ps, p)

            Puo = hist(Vuo, P)
            p = bar(Puo; label = "", title,
                xlabel = "uo", ylabel = "proj"
            )
            uo_av = sum(P .* Vuo)
            vline!(p, [uo_av]; label = "", color = :red)
            push!(ps, p)
        end

        layout = (2, 3)
        PltU.sfig(ps, 
            plotsdir(simid, "dyn_maxent_marginals", simparams, (;MEmode), ".png");
            layout
        )
    end

end
