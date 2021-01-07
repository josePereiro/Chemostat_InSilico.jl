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
function plot_res(M::SimModel, TS::ResTS; f = (x) -> x, marginf = 0.2)

    p1 = plot_ts_concs(TS; f, marginf)    
    p2 = plot_ts_D(TS; f, marginf)    
    p3 = plot_ts_X(TS; f, marginf)    

    p4 = plot_politope(M)
    
    p = plot([p1, p2, p3, p4]...;
        size = [800, 700], layout = 4)
end

## ----------------------------------------------------------------------------
function plot_ts_concs!(p, TS; f = (x) -> x, marginf = 0.2, pkwargs...)
    plot!(p, xlabel = "time", ylabel = "conc")
    if !(isempty(TS.sg_ts) || isempty(TS.sl_ts))
        ylim = lims(marginf, TS.sg_ts, TS.sl_ts)
        plot!(p, f.(TS.sg_ts); ylim, label = "sg", lw = 3)
        plot!(p, f.(TS.sl_ts); ylim, label = "sl", lw = 3)
    end  
    plot!(p; pkwargs...)  
end
plot_ts_concs(TS; f = (x) -> x, marginf = 0.2, pkwargs...) = 
    plot_ts_concs!(plot(), TS; f, marginf, pkwargs...)    

## ----------------------------------------------------------------------------
function plot_ts_D!(p, TS; f = (x) -> x, marginf = 0.2, pkwargs...)
    plot!(p; xlabel = "time", ylabel = "X")
    if !isempty(TS.X_ts)
        ylim = lims(marginf, TS.X_ts)
        plot!(p, f.(TS.X_ts); ylim, label = "X", lw = 3)
    end    
    plot!(p; pkwargs...)  
end
plot_ts_D(TS; f = (x) -> x, marginf = 0.2, pkwargs...) =  
    plot_ts_D!(plot(), TS; f, marginf, pkwargs...)    

## ----------------------------------------------------------------------------
function plot_ts_X!(p, TS; f = (x) -> x, marginf = 0.2, pkwargs...)
    plot!(p; xlabel = "time", ylabel = "D")
    if !isempty(TS.D_ts)
        ylim = lims(marginf, TS.D_ts)
        plot!(p, f.(TS.D_ts); ylim, label = "D", lw = 3)
    end
    plot!(p; pkwargs...)  
end
plot_ts_X(TS; f = (x) -> x, marginf = 0.2, pkwargs...) =  
    plot_ts_X!(plot(), TS; f, marginf, pkwargs...)    

## ----------------------------------------------------------------------------
function plot_marginals(marginals, rxns = keys(marginals))
    ps = []
    params = (;xlabel = "flx", ylabel = "prob", yaxis = nothing, grid = false, 
        titlefont = 10, xaxisfont = 10)
    for rxn in rxns
        marginal = marginals[rxn]
        # Plots
        try
            p = plot(;title = rxn, params...)
            plot!(p, marginal; label = "", color = :red, alpha = 0.9, lw = 3, ylim = [0.0, Inf])

            av = av(marginal)
            std = sqrt(va(marginal))

            vline!(p, [av]; label = "", alpha = 0.8, lw = 2, ls = :dash, color = :black)
            vline!(p, [av - std]; label = "", alpha = 0.6, lw = 2, ls = :dot, color = :black)
            vline!(p, [av + std]; label = "", alpha = 0.6, lw = 2, ls = :dot, color = :black)
            push!(ps, p)
        catch err 
            @error string("Doing $rxn\n", err_str(err))
        end
    end
    plot(ps...; layout = length(ps))
end

## ----------------------------------------------------------------------------
function plot_politope(M::SimModel; 
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