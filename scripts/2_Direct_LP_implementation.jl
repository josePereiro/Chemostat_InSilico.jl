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
sim_name(name; sim_params...) = string(InLP.mysavename(name, ""; sim_params...), "__", now())
fig_dir(sname) = joinpath(InLP.CACHE_DIR, sname, "figures")
dat_dir(sname) = joinpath(InLP.CACHE_DIR, sname, "dat")
function make_dirs(sname)
    fdir = fig_dir(sname)
    ddir = dat_dir(sname)
    mkpath.([fdir, ddir])
end

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

function push_plot_save(M, TS, it; 
        sname, sim_params,
        push_frec = 10,
        savefig_frec = 100,
        savedat_frec = 1000,
        force = false
    )

    if it % push_frec == 0 || force
        push!(TS, M)        
    end

    # save fig
    save_fig = it == 1 || rem(it, savefig_frec) == 0 || it == M.niters || force
    if save_fig
        fname = InLP.mysavename("fig", "png"; it, sim_params...)
        fpath = joinpath(fig_dir(sname), fname)
        p = InLP.plot_res(M, TS)
        savefig(p, fpath)

    end

    # save data
    save_dat = rem(it, savedat_frec) == 0 || it == M.niters || force
    if save_dat
        fname = InLP.mysavename("tot_dat", "jls"; sim_params...)
        fpath = joinpath(dat_dir(sname), fname)
        serialize(fpath, (;TS, M))
    end
end

## ---------------------------------------------------------
# Simulation 
sname = sim_name("LP_Implementation")
DAT = UJL.DictTree() # To Store relevant information
@time let
    # setup
    
    make_dirs(sname)

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

    @info "Starting simulation" sname

    # cache
    InLP.vgvatp_cache(M0)

    # simulation params
    stst_th = 0.05
    stst_window = 250
    check_stst_frec = 1000
    savefig_frec = 1000
    savedat_frec = 1000
    push_frec = 10
    dead_th = DAT[:dead_th] = 1e-2

    # Test
    ϵs = [0.02:0.02:0.1; 0.2:0.1:1.0]
    Ds = [0.0; 0.001:0.001:0.015; 10.0.^-(1.6:0.1:2.2)] |> unique |> sort
    # Vls = [0.0, 0.1]
    Vls = [0.0]
    # τs = [0.0, 0.0022]
    τs = [0.0] |> sort
    
    c = 1
    N = prod(length.([ϵs, Ds, Vls, τs]))
    @info "Computing $(N) iterations"

    for D in Ds, ϵ in ϵs, Vl in Vls, τ in τs

        M = deepcopy(M0)
        M.D, M.ϵ, M.Vl, M.τ = D, ϵ, Vl, τ
        TS = InLP.ResTS()

        @info "Doing $c/$N ..." D ϵ Vl τ now()

        function on_iter(it, M)

            sim_params = (;M.D, M.ϵ, M.Δt, M.σ, M.τ, M.cg, M.Vl)

            # stead state
            stst = false
            if rem(it, check_stst_frec) == 0
                stst = check_stst(TS; w = stst_window, th = stst_th)
            end
            
            # cells die
            dead = M.X < dead_th

            # finish
            finish = dead || stst || it == M.niters

            # output
            push_plot_save(M, TS, it; 
                sname, sim_params,
                push_frec, savefig_frec, savedat_frec,
                force = finish
            )

            return finish
        end
        InLP.run_simulation!(M; on_iter, verbose = true, force = false)
    end
end

## ----------------------------------------------------------------------------
# Collecting Bundle
let
    ddir = dat_dir(sname) 
    @assert isdir(ddir)

    dfiles = readdir(ddir)
    N  = length(dfiles)
    for (i, file) in dfiles |> enumerate
        TS, M = deserialize(joinpath(ddir, file))
        @info "Doing $i/$N ... " M.D M.ϵ file; println()
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
UJL.save_data(InCh.CH3_DAT_BUNDLE_FILE, DAT)
