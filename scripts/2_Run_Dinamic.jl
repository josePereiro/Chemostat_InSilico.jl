import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    using Base.Threads
    using Dates
    using Statistics
    using InteractiveUtils
    using Random

end

## -----------------------------------------------------------------------------------------------
# Tools
const WLOCK = ReentrantLock()
const DATA_DIR = InCh.DYN_DATA_DIR
const DATA_FILE_PREFFIX = "dyn_dat"
dat_file(;sim_params...) = joinpath(DATA_DIR, 
    InLP.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...))

function check_stst(ts; stst_window, stst_th)
    @views for sym in [:X_ts, :sl_ts, :sg_ts, :D_ts]
        s = getfield(ts, sym)
        length(s) <= stst_window && return false
        vs = s[end - stst_window:end]
        m = abs(mean(vs))
        m = m == 0 ? 1.0 : m
        s = std(vs)
        s/m > stst_th && return false
    end
    return true
end

# ## ----------------------------------------------------------------------------
# let
#     M0 = InLP.SimModel(;
#             δvatp = 2, 
#             δvg = 3, 
#             niters = 20,
#             X0 = 1.3,
#             sg0 = 15.0,
#             sl0 = 0.0,
#             Δt = 0.5,
#         )

#     # on_iter(it, M) = (sleep(1.0); false)
#     on_iter(it, M) = false
#     InLP.run_simulation!(M0; on_iter, verbose_frec = 10)
# end;

# ## ----------------------------------------------------------------------------
# ## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
# Simulation 
INDEX = UJL.DictTree() # To Store relevant information
@time let
    
    # setup
    mkpath(DATA_DIR)

    # base model
    M0 = INDEX[:M0] = InLP.SimModel(;
            δvatp = 2, 
            δvg = 3, 
            niters = Int(5e4),
            X0 = 0.3,
            sg0 = 15.0,
            sl0 = 0.0,
            Δt = 0.5,
        )

    @info "Starting simulation" nthreads() now()
    println()

    # cache
    LP_cache = InLP.vgvatp_cache(M0)

    # simulation params
    stst_th = 0.05
    stst_window = 250
    check_stst_frec = 1000
    savedat_frec = 1000
    info_frec = 100
    push_frec = 10
    death_th = INDEX[:death_th] = 1e-2

    # Params
    # Vls = [0.0, 0.1]
    Vls = INDEX[:Vls] = [0.0]
    # Ds = INDEX[:Ds]= [0.003:0.001:0.045;]
    Ds = INDEX[:Ds] = [0.003:0.005:0.045;]
    ϵs = INDEX[:ϵs] = [0.01, 0.1, 0.3, 0.5, 0.8, 1.0]
    # τs = [0.0, 0.0022]
    τs = INDEX[:τs] = [0.0]
    
    # 
    params = Iterators.product(Vls, Ds, ϵs, τs)
    N = length(params)
    @info("Computing $(N) iterations"); println() 
    gc = 0
    
    # feeding task
    Ch = Channel(1) do Ch_
        for (Vl, D, ϵ, τ) in params
            put!(Ch_, (Vl, D, ϵ, τ))
        end 
    end

    @threads for thid in 1:nthreads()
        for (Vl, D, ϵ, τ) in Ch

            # This must identify the iteration
            sim_params = (;D, ϵ, τ, Vl, M0.Δt, M0.σ, M0.cg)

            # Saying hello
            c = nothing
            lock(WLOCK) do
                gc += 1; c = gc
                @info("Starting $c/$N ...", 
                    (Vl, D, ϵ, τ), 
                    now(), thid
                ); println()
            end
            
            # Setup or load cache
            cfile = dat_file(;sim_params...)
            if isfile(cfile)
                status, TS, M = deserialize(cfile)
                tslen = length(TS.X_ts)
                lock(WLOCK) do
                    @info("Cache loaded $c/$N ... ", 
                        M.X, M.sl, M.sg,
                        tslen, status, 
                        (Vl, D, ϵ, τ), 
                        now(), thid
                    ); println()
                end
            else
                # setup model
                M = deepcopy(M0)
                M.D, M.ϵ, M.Vl, M.τ = D, ϵ, Vl, τ
                TS = InLP.ResTS()
                tslen = length(TS.X_ts)
                status = :running
            end
            
            # Test
            status == :finished && (status = :running) # Compat

            function on_iter(it, _)

                # update tslen
                tslen = length(TS.X_ts)
                
                # check steady state
                stst = rem(it, check_stst_frec) == 0 && 
                    check_stst(TS; stst_window, stst_th)
                stst && (status = :stst)
                
                # check cells die
                death = M.X < death_th
                death && (status = :death)

                # simulation ends
                sim_ends = it == M.niters
                sim_ends && (status = :sim_ends)

                # finish
                finished = death || stst || sim_ends

                push_state = rem(it, push_frec) == 0 || finished
                push_state && push!(TS, M)        

                # save data
                save_dat = rem(it, savedat_frec) == 0 || finished
                save_dat && serialize(cfile, (;status, TS, M))
                save_dat && touch(cfile) # plotting trigger
                    
                # info and lock dat
                show_info = rem(it, info_frec) == 0
                if show_info
                    lock(WLOCK) do
                        @info("Doing $c/$N ... ", 
                            it, M.X, M.sl, M.sg, 
                            tslen, status, 
                            (Vl, D, ϵ, τ),
                            now(), thid
                        ); println()
                    end
                end

                return finished
            end
            
            if status == :running 
                InLP.run_simulation!(M; 
                    on_iter, LP_cache, 
                    verbose = false
                )
            end
            
            lock(WLOCK) do
                @info("Finished $c/$N ... ", 
                    M.X, M.sl, M.sg, 
                    tslen, status, 
                    (Vl, D, ϵ, τ), 
                    now(), thid
                ); println()
                INDEX[:DFILE, Vl, D, ϵ, τ] = cfile
                GC.gc()
            end
        end # for (Vl, D, ϵ, τ)
    end # @threads for thid
end

## -----------------------------------------------------------------------------------------------
# SAVING
UJL.save_data(InCh.DYN_DATA_INDEX_FILE, INDEX)
