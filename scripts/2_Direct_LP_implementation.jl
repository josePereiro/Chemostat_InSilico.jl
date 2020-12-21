import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

import Chemostat_InSilico
const InLP = Chemostat_InSilico.LP_Implement
const InU = Chemostat_InSilico.Utilities

using ProgressMeter
using Plots
using Serialization
using Base.Threads
using Dates
using Base.Threads
using BenchmarkTools
using Statistics

## ---------------------------------------------------------
# Tools
sim_name(name) = string(name, "__", now())
fig_dir(sname) = joinpath(InLP.CH3_FIGURES_DIR, sname)
dat_dir(sname) = joinpath(InLP.CH3_DATA_DIR, sname)
function make_dirs(sname)
    fdir = fig_dir(sname)
    ddir = dat_dir(sname)
    mkpath.([fdir, ddir])
end

function check_stst(ts; w = 1000, th = 0.05)
    stst = true
    @views for sym in [:X_ts, :sl_ts, :sg_ts, :D_ts]
        s = getfield(ts, sym)
        length(s) <= w && return false
        vs = s[end - w:end]
        m = abs(mean(vs))
        m = m == 0 ? 1.0 : m
        s = std(vs)
        stst = stst && s/m < th
    end
    return stst
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
# Common model
M0 = InLP.SimModel(;
            θvatp = 2, 
            θvg = 3, 
            niters = Int(5e4),
            X0 = 1.5,
            sg0 = 15.0,
            sl0 = 2.0,
            damp = 0.0,
            Δt = 0.5,
            ϵ = 0.5
        )
## ---------------------------------------------------------
let

    # setup
    sname = sim_name("Ch3_")
    make_dirs(sname)

    @info "Doing" sname

    # cache
    InLP.vgvatp_cache(M0)

    # simulation params
    stst_th = 0.05
    stst_window = 250
    check_stst_frec = 1000
    savefig_frec = 500
    savedat_frec = 1000
    push_frec = 10
    dead_th = 1e-2

    Ds = collect(0.001:0.001:0.1) |> sort
    # Ds = [Ds; reverse(Ds)] |> unique
    for D in Ds
        M = deepcopy(M0)
        TS = InLP.ResTS()
        M.D = D

        function on_iter(it, M)

            sim_params = (;M.D, M.ϵ, M.Δt)

            # stead state
            stst = false
            if rem(it, check_stst_frec) == 0
                stst = check_stst(TS; w = stst_window, th = stst_th)
            end
            
            # save stst
            # if stst
            #     fname = InLP.mysavename("stst_dat", "jls"; it, sim_params...)
            #     fpath = joinpath(dat_dir(sname), fname)
            #     serialize(fpath, (;it, M))
            # end

            # cells die
            dead = M.X < dead_th

            # change D
            # Ds_empty = isempty(Ds) 
            # stst && !Ds_empty && (M.D = popfirst!(Ds))

            # finish
            # finish = dead || (stst && Ds_empty)
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
