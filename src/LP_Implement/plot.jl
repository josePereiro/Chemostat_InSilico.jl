function lims(marginf, s...)
    m = minimum(minimum.(s))
    M = maximum(maximum.(s))
    gamma = abs(m - M)
    gamma = gamma == 0 ? M : gamma
    gamma = gamma == 0 ? marginf : gamma
    margin = gamma * marginf
    (m - margin, M + margin)
end

## ----------------------------------------------------------------------------
function plot_res(M::SimModel, ts::ResTS; f = (x) -> x, marginf = 0.2)

    p1 = plot_ts_concs(ts; f, marginf)    
    p2 = plot_ts_D(ts; f, marginf)    
    p3 = plot_ts_X(ts; f, marginf)    

    p4 = plot_politope(M)
    
    p = plot([p1, p2, p3, p4]...;
        size = [800, 700], layout = 4)
end

## ----------------------------------------------------------------------------
function plot_ts_concs!(p, ts; f = (x) -> x, marginf = 0.2, pkwargs...)
    plot!(p, xlabel = "time", ylabel = "conc")
    if !(isempty(ts.sg_ts) || isempty(ts.sl_ts))
        ylim = lims(marginf, ts.sg_ts, ts.sl_ts)
        plot!(p, f.(ts.sg_ts); ylim, label = "sg", lw = 3)
        plot!(p, f.(ts.sl_ts); ylim, label = "sl", lw = 3)
    end  
    plot!(p; pkwargs...)  
end
plot_ts_concs(ts; f = (x) -> x, marginf = 0.2, pkwargs...) = 
    plot_ts_concs!(plot(), ts; f, marginf, pkwargs...)    

## ----------------------------------------------------------------------------
function plot_ts_D!(p, ts; f = (x) -> x, marginf = 0.2, pkwargs...)
    plot!(p; xlabel = "time", ylabel = "X")
    if !isempty(ts.X_ts)
        ylim = lims(marginf, ts.X_ts)
        plot!(p, f.(ts.X_ts); ylim, label = "X", lw = 3)
    end    
    plot!(p; pkwargs...)  
end
plot_ts_D(ts; f = (x) -> x, marginf = 0.2, pkwargs...) =  
    plot_ts_D!(plot(), ts; f, marginf, pkwargs...)    

## ----------------------------------------------------------------------------
function plot_ts_X!(p, ts; f = (x) -> x, marginf = 0.2, pkwargs...)
    plot!(p; xlabel = "time", ylabel = "D")
    if !isempty(ts.D_ts)
        ylim = lims(marginf, ts.D_ts)
        plot!(p, f.(ts.D_ts); ylim, label = "D", lw = 3)
    end
    plot!(p; pkwargs...)  
end
plot_ts_X(ts; f = (x) -> x, marginf = 0.2, pkwargs...) =  
    plot_ts_X!(plot(), ts; f, marginf, pkwargs...)    

## ----------------------------------------------------------------------------
function plot_marginals(marginals, rxns = keys(marginals))
    ps = []
    params = (;xlabel = "flx", ylabel = "prob", yaxis = nothing, grid = false)
    for rxn in rxns
        marginal = marginals[rxn]
        # Plots
        try
            p = plot(;title = rxn, params...)
            plot!(p, marginal; label = "", color = :red, alpha = 0.9, lw = 3, ylim = [0.0, Inf])

            av = marginal_av(marginal)
            std = sqrt(marginal_va(marginal))

            vline!(p, [av]; label = "", alpha = 0.7, lw = 2, ls = :dash, color = :black)
            vline!(p, [av - std]; label = "", alpha = 0.3, lw = 2, ls = :dot, color = :black)
            vline!(p, [av + std]; label = "", alpha = 0.3, lw = 2, ls = :dot, color = :black)
            push!(ps, p)
        catch err 
            @error string("Doing $rxn\n", err_str(err))
        end
    end
    plot(ps...; layout = length(ps))
end

## ----------------------------------------------------------------------------
function plot_politope(M::SimModel; 
        D = 250.0,
        hits_count = 3000,
        maxiters = 100 * hits_count,
    )
    
    vatp_range, vg_ranges = vatpvg_ranges(M)
    vgL, vgU = minimum(first.(vg_ranges)), maximum(last.(vg_ranges))
    Δvg = step(first(vg_ranges))
    vg_range = vgL:Δvg:vgU
    mX, MX = lXgamma(M)
    vatpvgN = sum(length.(vg_ranges))

    p = plot(xlabel = "vatp", ylabel = "vg")
    hits = 0
    iters = 0
    while vatpvgN > 1 && hits < hits_count && iters < maxiters
        iters += 1
        vatp = rand(vatp_range)
        vg = rand(vg_range)
        !haskey(M.Xb, vatp) && continue
        !haskey(M.Xb[vatp], vg) && continue

        lX = M.Xb[vatp][vg]
        color = :black
        ms = 10.0 * lX/MX
        scatter!(p, [vatp], [vg]; color, label = "", ms, alpha = 0.2)
        hits += 1
    end
    p
end