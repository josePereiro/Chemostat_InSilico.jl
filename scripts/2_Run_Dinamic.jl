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
    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    using Base.Threads
    using Dates
    using Statistics
    using InteractiveUtils

end

## ---------------------------------------------------------
# Tools
PROG_FIG_DIR = joinpath(InLP.DYN_FIGURES_DIR, "progress")
CACHE_DIR = InLP.DYN_CACHE_DIR
dat_file(;sim_params...) = joinpath(CACHE_DIR, 
    InLP.mysavename("dyn_dat", "jls"; sim_params...))
fig_file(;sim_params...) = joinpath(PROG_FIG_DIR, 
    InLP.mysavename("fig", "png"; sim_params...))

make_dirs() = mkpath.([PROG_FIG_DIR, CACHE_DIR])

function check_stst(ts; w = 1000, th = 0.05)
    @views for sym in [:X_ts, :sl_ts, :sg_ts, :D_ts]
        s = getfield(ts, sym)
        length(s) <= w && return false
        vs = s[end - w:end]
        m = abs(mean(vs))
        m = m == 0 ? 1.0 : m
        s = std(vs)
        s/m > th && return false
    end
    return true
end

## ---------------------------------------------------------
# Simulation 
DAT = UJL.DictTree() # To Store relevant information
@time let
    # setup
    
    make_dirs()

    # base model
    M0 = InLP.SimModel(;
            δvatp = 2, 
            δvg = 3, 
            niters = Int(5e4),
            X0 = 1.5,
            sg0 = 15.0,
            sl0 = 0.0,
            Δt = 0.5,
        )

    @info "Starting simulation" now()

    # cache
    InLP.vgvatp_cache(M0)

    # simulation params
    stst_th = 0.05
    stst_window = 250
    check_stst_frec = 1000
    savefig_frec = 1000
    savedat_frec = 1000
    push_frec = 10
    death_th = DAT[:death_th] = 1e-2

    # Params
    ϵs = [0.01; 0.1:0.1:1.0]
    Ds = [0.001:0.001:0.017;]
    # Vls = [0.0, 0.1]
    Vls = [0.0]
    # τs = [0.0, 0.0022]
    τs = [0.0] 
    
    c = 1
    N = prod(length.([ϵs, Ds, Vls, τs]))
    @info "Computing $(N) iterations"
    cfiles = Set([])

    for D in Ds, ϵ in ϵs, Vl in Vls, τ in τs

        # Prepare simulation model
        M = deepcopy(M0)
        M.D, M.ϵ, M.Vl, M.τ = D, ϵ, Vl, τ
        TS = InLP.ResTS()
        sim_params = (;M.D, M.ϵ, M.Δt, M.σ, M.τ, M.cg, M.Vl)

        @info "Doing $c/$N ..." D ϵ Vl τ now()s
        
        function on_iter(it, _)

            # check cache
            cfile = dat_file(;sim_params...)
            push!(cfiles, cfile)
            cached = it == 1 && isfile(cfile)
            cached && (TS, M = deserialize(cfile))
            
            # check steady state
            stst = rem(it, check_stst_frec) == 0 && 
                check_stst(TS; w = stst_window, th = stst_th)
            
            # check cells die
            death = M.X < death_th

            # finish
            finish = death || stst || it == M.niters

            push_state = it % push_frec == 0 || finish
            push_state && push!(TS, M)        

            # save fig
            save_fig = it == 1 || rem(it, savefig_frec) == 0 || 
                it == M.niters || finish || !cached
            save_fig && savefig(InLP.plot_res(M, TS), fig_file(;it, sim_params...))

            # save data
            save_dat = rem(it, savedat_frec) == 0 || it == M.niters || 
                finish || !cached
            save_dat && serialize(cfile, (;TS, M))

            return finish
        end
        InLP.run_simulation!(M; on_iter, verbose = true, force = false)

        c += 1
    end

    ## ----------------------------------------------------------------------------
    # Collecting and Bundle    
    for (i, cfile) in enumerate(cfiles)
        TS, M = deserialize(cfile)
        @info "Doing $i/$N ... " M.D M.ϵ cfile; println()
        DAT[:TS, M.Vl, M.D, M.ϵ] = TS
        DAT[:M, M.Vl, M.D, M.ϵ] = M

        push!(get!(DAT, [], :Ds), M.D)
        push!(get!(DAT, [], :Vls), M.Vl)
        push!(get!(DAT, [], :ϵs), M.ϵ)
    end
    sort!(unique!(DAT[:Ds]))
    sort!(unique!(DAT[:Vls]))
    sort!(unique!(DAT[:ϵs]))
end
varinfo(Main, r"DAT")

## ----------------------------------------------------------------------------
# SAVING
UJL.save_data(InCh.DYN_DATA_BUNDLE_FILE, DAT)
