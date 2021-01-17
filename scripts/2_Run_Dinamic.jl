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

## ---------------------------------------------------------
# Tools
const WLOCK = ReentrantLock()
PROG_FIG_DIR = joinpath(InCh.DYN_FIGURES_DIR, "progress")
DATA_DIR = InCh.DYN_DATA_DIR
DATA_FILE_PREFFIX = "dyn_dat"
dat_file(;sim_params...) = joinpath(DATA_DIR, 
    InLP.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...))
fig_file(;sim_params...) = joinpath(PROG_FIG_DIR, 
    InLP.mysavename("fig", "png"; sim_params...))

make_dirs() = mkpath.([PROG_FIG_DIR, DATA_DIR])

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

## ---------------------------------------------------------
# Plot progress coroutine
const MTIMES = Dict()
function plot_progress(mtimes = MTIMES; force = false)
    for file in readdir(DATA_DIR)
        !startswith(file, DATA_FILE_PREFFIX) && continue
        
        cfile = joinpath(DATA_DIR, file)

        currt = mtime(cfile)
        lastt = get!(mtimes, file, currt)

        # save fig
        if lastt != currt || force
            status, TS, M = deserialize(cfile)
            p = InLP.plot_res(M, TS)
            fname = replace(file, DATA_FILE_PREFFIX => "fig")
            fname, _ = splitext(fname)
            fname = string(fname, ".png")

            lock(WLOCK) do
                @info "Plotting progress" fname now() threadid()
                println()
            end
            
            ffile = joinpath(PROG_FIG_DIR, fname)
            savefig(p, ffile)
            mtimes[file] = currt

        end
    end
end
@async let
    sleep(120) # Do not interfere at start
    while true
        sleep(rand(30:60))
        plot_progress()
    end
end

## ---------------------------------------------------------
# Simulation 
INDEX = UJL.DictTree() # To Store relevant information
@time let
    
    # setup
    make_dirs()

    # base model
    M0 = INDEX[:M0] = InLP.SimModel(;
            δvatp = 2, 
            δvg = 3, 
            niters = Int(5e4),
            X0 = 1.5,
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
    Ds = INDEX[:Ds]= [0.003:0.001:0.045;]
    ϵs = INDEX[:ϵs] = [0.01, 0.1, 0.5, 1.0]
    # τs = [0.0, 0.0022]
    τs = INDEX[:τs] = [0.0] 
    
    params = Iterators.product(Vls, Ds, ϵs, τs) |> collect |> shuffle! # distribute equally between threads
    N = length(params)
    @info "Computing $(N) iterations"
    println() 

    c = 0
    @threads for (Vl, D, ϵ, τ) in params

        # This must identify the iteration
        sim_params = (;D, ϵ, τ, Vl, M0.Δt, M0.σ, M0.cg)

        # Say hello
        lock(WLOCK) do
            c += 1
            @info "Starting $c/$N ..." Vl, D, ϵ, τ now() threadid()
            println()
        end
        
        # Setup or load cache
        cfile = dat_file(;sim_params...)
        if isfile(cfile)
            status, TS, M = deserialize(cfile)
            tslen = length(TS.X_ts)
            lock(WLOCK) do
                @info "Cache loaded $c/$N ... " M.X M.sl M.sg tslen status Vl, D, ϵ, τ now() threadid()
                println()
            end
        else
            M = deepcopy(M0)
            M.D, M.ϵ, M.Vl, M.τ = D, ϵ, Vl, τ
            TS = InLP.ResTS()
            tslen = length(TS.X_ts)
            status = :running
        end
        
        # Test
        status == :finished && (status = :running) # Compat

        if status == :running
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
                        @info "Doing $c/$N ... " it M.X M.sl M.sg tslen status Vl, D, ϵ, τ now() threadid()
                        println()
                    end
                end

                return finished
            end
            InLP.run_simulation!(M; on_iter, LP_cache, verbose = false, force = false)
        end
        
        lock(WLOCK) do
            @info "Finished $c/$N ... " M.X M.sl M.sg tslen status Vl, D, ϵ, τ now() threadid()
            INDEX[:DFILE, Vl, D, ϵ, τ] = cfile
            println()
            GC.gc()
        end
    end
    plot_progress(;force = true)
end

## ----------------------------------------------------------------------------
# SAVING
UJL.save_data(InCh.DYN_DATA_INDEX_FILE, INDEX)