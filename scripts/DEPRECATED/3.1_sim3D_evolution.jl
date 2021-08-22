## ------------------------------------------------------
# evolution gif
let
    # Test
    # SimD3 [D=3.29e-01 cg=1.50e+01 ϵ=1.00e+00].status.jls
    # simid = "SimD3"
    # Ds = [3.29e-01]
    # cgs = [1.50e+01] 
    # ϵs = [1.00e+00]
    
    nframes = 30
    plk = ReentrantLock()
    iter = enumerate(Iterators.product(Ds, ϵs, cgs))
    @threads for (simi, (D, ϵ, cg)) in collect(iter)

        thid = threadid()
       
        status = get_status(simid, (;D, ϵ, cg))
        (status == UNDONE_SIM_STATUS) && continue
        
        (ts, Pzts)    = collect_ts(:Pzts, simid, (;D, ϵ, cg))
        (ts, Pugts)   = collect_ts(:Pugts, simid, (;D, ϵ, cg))
        (ts, Puots)   = collect_ts(:Puots, simid, (;D, ϵ, cg))
        (ts, Xts)     = collect_ts(:Xts, simid, (;D, ϵ, cg))
        (ts, z_avts)  = collect_ts(:z_avts, simid, (;D, ϵ, cg))
        (ts, ug_avts) = collect_ts(:ug_avts, simid, (;D, ϵ, cg))
        (ts, uo_avts) = collect_ts(:uo_avts, simid, (;D, ϵ, cg))
        (ts, cgD_Xts) = collect_ts(:cgD_Xts, simid, (;D, ϵ, cg))
        (ts, sgts)    = collect_ts(:sgts, simid, (;D, ϵ, cg))

        @info("At", D, ϵ, cg, status, thid, length(ts))
        isempty(ts) && continue

        prms = (;label = "", xrotation = 45)
        barprms = (;yticks = false)
        plotsprms = (;lw = 3, ylim = (0, Inf))
        
        idxs = firstindex(ts):max(1, div(length(ts), nframes)):lastindex(ts)
        gifs_ = lock(plk) do
            return map(enumerate(idxs)) do (i, idx)

                title = string(status, " at ", round(ts[idx]; sigdigits = 3))
            
                # hist
                pz_hist = bar(Pzts[idx]; 
                    xlabel = "z", ylabel = "proj", title, 
                    barprms...,
                    prms...
                )
                pug_hist = bar(Pugts[idx]; 
                    xlabel = "ug", ylabel = "proj", title, 
                    barprms...,
                    prms...
                )
                puo_hist = bar(Puots[idx]; 
                    xlabel = "uo", ylabel = "proj", title, 
                    barprms...,
                    prms...
                )

                # time series
                ts_idxs = idxs[begin:i]
                pX_ts = plot(ts[ts_idxs], Xts[ts_idxs]; 
                    xlabel = "time", ylabel = "X", title, 
                    plotsprms..., 
                    prms...
                )
                
                pz_avts = plot(ts[ts_idxs], z_avts[ts_idxs]; 
                    xlabel = "time", ylabel = "z_av & D", title, 
                    plotsprms..., 
                    prms...
                )
                plot!(pz_avts, ts[ts_idxs], fill(D, length(ts_idxs)); 
                    color = :red, ls = :dash, 
                    plotsprms..., 
                    prms...
                )

                p_sgts = plot(ts[ts_idxs], sgts[ts_idxs]; 
                    xlabel = "time", ylabel = "sgts & cg", title, 
                    plotsprms..., 
                    prms...
                )
                plot!(p_sgts, ts[ts_idxs], fill(cg, length(ts_idxs)); 
                    color = :red, ls = :dash, 
                    plotsprms..., 
                    prms...
                )
                
                pug_avts = plot(ts[ts_idxs], ug_avts[ts_idxs]; 
                    xlabel = "time", ylabel = "ug_av & cgD/X", title, 
                    plotsprms..., 
                    prms...
                )
                plot!(pug_avts, ts[ts_idxs], cgD_Xts[ts_idxs]; 
                    color = :red, ls = :dash, 
                    plotsprms..., 
                    prms...
                )

                puo_avts = plot(ts[ts_idxs], uo_avts[ts_idxs]; 
                    xlabel = "time", ylabel = "uo_av", title, 
                    plotsprms..., 
                    prms...
                )

                all_plots = [
                    pz_hist, pug_hist, puo_hist, 
                    pz_avts, pug_avts, puo_avts, pX_ts, p_sgts
                ]
                PltU.sfig(all_plots, tempname(), ".png")
            end
        end  # lock

        # last frame
        cp(
            last(gifs_), 
            plotsdir(simid, "dyn_evolution", (;D, ϵ, cg), ".png");
            force = true
        )
        
        # make gif (slow part)
        PltU.save_gif(gifs_, 
            plotsdir(simid, "dyn_evolution", (;D, ϵ, cg), ".gif")
        )

        rm.(gifs_; force = true)

    end # @threads for
end