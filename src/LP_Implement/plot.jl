function lims(marginf, s...)
    m = minimum(minimum.(s))
    M = maximum(maximum.(s))
    gamma = abs(m - M)
    gamma = gamma == 0 ? M : gamma
    gamma = gamma == 0 ? marginf : gamma
    margin = gamma * marginf
    (m - margin, M + margin)
end

## ---------------------------------------------------------
function plot_res(M::SimModel, ts::ResTS; f = (x) -> x, marginf = 0.2)
    
    p1 = plot(xlabel = "time", ylabel = "conc")
    if !(isempty(ts.sg_ts) || isempty(ts.sl_ts))
        ylim = lims(marginf, ts.sg_ts, ts.sl_ts)
        plot!(p1, f.(ts.sg_ts); ylim, label = "sg", lw = 3)
        plot!(p1, f.(ts.sl_ts); ylim, label = "sl", lw = 3)
    end    
    
    p2 = plot(xlabel = "time", ylabel = "X")
    if !isempty(ts.X_ts)
        ylim = lims(marginf, ts.X_ts)
        plot!(p2, f.(ts.X_ts); ylim, label = "X", lw = 3)
    end    
    
    p3 = plot(xlabel = "time", ylabel = "D")
    if !isempty(ts.D_ts)
        ylim = lims(marginf, ts.D_ts)
        plot!(p3, f.(ts.D_ts); ylim, label = "D", lw = 3)
    end

    p4 = plot_politope(M)
    
    p = plot([p1, p2, p3, p4]...;
        size = [800, 700], layout = 4)
end

function plot_politope(M::SimModel; 
        D = 250.0,
        N = 3000
    )
    
    vatp_range, vg_ranges = vatpvg_ranges(M)
    vgL, vgU = minimum(first.(vg_ranges)), maximum(last.(vg_ranges))
    Δvg = step(first(vg_ranges))
    vg_range = vgL:Δvg:vgU
    mX, MX = lXgamma(M)

    p = plot(xlabel = "vatp", ylabel = "vg")
    c = 0
    while c < N
        vatp = rand(vatp_range)
        vg = rand(vg_range)
        !haskey(M.Xb, vatp) && continue
        !haskey(M.Xb[vatp], vg) && continue

        lX = M.Xb[vatp][vg]
        color = :black
        ms = 10.0 * lX/MX
        scatter!(p, [vatp], [vg]; color, label = "", ms, alpha = 0.2)
        c += 1
    end
    p
end