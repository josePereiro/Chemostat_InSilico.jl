using ProjAssistant
@quickactivate

## ------------------------------------------------------
@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        maxent, normalizeP!, load_simdat, get_status,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS, collect_ts
        
    using Base.Threads
    using ExtractMacro

end

## ----------------------------------------------------------------------------
# ME_AVZ_EXPECTED_AVUG_BOUNDED
let
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    # batch = (; push_frec, Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts)
    simid = :SimD3
    params = lglob(Dyn, simid, :params, :finite_cg)
    @extract params: Ds ϵs cg

    _simparams_pool = Iterators.product(Ds, ϵs, cg)
    tot_sims = length(_simparams_pool)
    Ch = Channel() do Ch_
        for (simi, (D, ϵ, cg)) in enumerate(_simparams_pool)
            put!(Ch_, (simi, (D, ϵ, cg)))
        end
    end

    @threads for _ in 1:nthreads()
        thid = threadid()

        for (simi, (D, ϵ, cg)) in Ch
            simparams = (;D, ϵ, cg)

            # ---------------------------------------------------------------
            # Welcome && status && cache
            status = Dyn.get_status(simid, (;D, ϵ, cg))
            datfile = Dyn.maxent_file(simid, MEmode, simparams)
            
            @info("At", simid, simparams, thid, status, relpath(datfile))
            isfile(datfile) && continue
            (status != STST_SIM_STATUS) && continue

            # ---------------------------------------------------------------
            # Sim globals
            S = Dyn.load_simdat(simid, simparams)
            X = S.X 
            cgD_X = S.cg * S.D / S.X
            V = S.V
            S = nothing; GC.gc()

            # ---------------------------------------------------------------
            # maxent
            me_out = maxent(V, D, cgD_X; 
                # Top loop
                top_maxiter = 100,

                # GD globals
                gd_maxiter = 1000,
                gdth = 1e-2,
                verbose = true,

                # stst
                stst_th = 1e-2,
                stst_w = 10, 

                # info
                further_info = (;
                    progress = string(simi, " of ", tot_sims, " (", div(simi * 100, tot_sims), "%)")
                )
            )

            sdat(Dyn, me_out, datfile)
        
        end # for (simi, (D, ϵ, cg))    

    end # for thid
end
