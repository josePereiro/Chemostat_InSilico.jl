using ProjAssistant
@quickactivate

## ------------------------------------------------------
@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        maxent, normalizeP!, load_simdat, get_status,
        load_batch, save_batch,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS, collect_ts,
        STST, is_steady
        

    using Plots
    import GR
    !isinteractive() && GR.inline("png")

    using Base.Threads
    using ProgressMeter
    using BenchmarkTools
    using ExtractMacro

end

## ------------------------------------------------------
# globlas
# batch = (; push_frec, Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts)
params = lglob(Dyn, :dyn, :params, :finite_cg)
@extract params: Ds ϵs cg simid

simparams = let
    iter = collect(Iterators.product(Ds, ϵs, cg))
    filter(iter) do tp
        D, ϵ, cg = tp
        status = Dyn.get_status(simid, (;D, ϵ, cg))
        (status != STST_SIM_STATUS) && return false
        return true
    end
end

## ----------------------------------------------------------------------------
# ME_AVZ_EXPECTED_AVUG_BOUNDED
let

    WLOCK = ReentrantLock()
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    sim_tot = length(simparams)
    @threads for (simi, (D, ϵ, cg)) in collect(enumerate(simparams))
        thid = threadid()

        simparams = (;D, ϵ, cg)

        # ---------------------------------------------------------------
        # Welcome && status
        @info("At", simid, simparams, thid)

        # ---------------------------------------------------------------
        datfile = Dyn.maxent_file(simid, MEmode, simparams)
        isfile(datfile) && continue

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
            stst_w = 10
        )

        sdat(Dyn, me_out, datfile)
        
    end # for (simi, (D, ϵ, cg)) 
end

## ----------------------------------------------------------------------------
let
    net3D = Dyn.ToyModel3D()
end