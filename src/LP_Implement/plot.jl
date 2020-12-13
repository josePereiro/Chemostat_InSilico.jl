## ---------------------------------------------------------
function plot_res(sg_ts, sl_ts, X_ts, D_ts; f = (x) -> x)
    p1 = plot(xlabel = "time", ylabel = "conc")
    plot!(p1, f.(sg_ts); label = "sg", lw = 3)
    plot!(p1, f.(sl_ts); label = "sl", lw = 3)
    
    p2 = plot(xlabel = "time", ylabel = "X")
    plot!(p2, f.(X_ts); label = "X", lw = 3)

    p3 = plot(xlabel = "time", ylabel = "D")
    plot!(p3, f.(D_ts); label = "D", lw = 3)
    
    p = plot([p1, p2, p3]...;
        size = [1000, 400], layout = grid(1, 3))
end
    
plot_res(ts::ResTS; f = (x) -> x) = plot_res(ts.sg_ts, ts.sl_ts, ts.X_ts, ts.D_ts; f)