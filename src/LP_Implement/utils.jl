## ---------------------------------------------------------
function plot_res(M::SimModel; f = (x) -> x)
    p1 = plot(xlabel = "time", ylabel = "conc")
    plot!(p1, f.(M.sg_ts); label = "sg", lw = 3)
    plot!(p1, f.(M.sl_ts); label = "sl", lw = 3)
    
    p2 = plot(xlabel = "time", ylabel = "X")
    plot!(p2, f.(M.X_ts); label = "X", lw = 3)
    
    p = plot([p1, p2]...)
end
    