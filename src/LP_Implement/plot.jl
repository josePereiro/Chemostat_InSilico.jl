## ----------------------------------------------------------------------------
minmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

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

    # time series
    p1 = plot_ts_concs(TS; f, marginf)    
    p2 = plot_ts_D(TS; f, marginf)    
    p3 = plot_ts_X(TS; f, marginf)    

    # polytope
    p4 = plot(xlabel = "vatp", ylabel = "vg")
    plot_polborder!(p4, M)
    plot_poldist!(p4, M; hits_count = 500,
        static_th = 0.05, rand_th = 1.0
    )
    
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

            av_ = av(marginal)
            std_ = sqrt(va(marginal))

            vline!(p, [av_]; label = "", alpha = 0.8, lw = 2, ls = :dash, color = :black)
            vline!(p, [av_ - std_]; label = "", alpha = 0.6, lw = 2, ls = :dot, color = :black)
            vline!(p, [av_ + std_]; label = "", alpha = 0.6, lw = 2, ls = :dot, color = :black)
            push!(ps, p)
        catch err 
            @error string("Doing $rxn\n", err_str(err))
        end
    end
    plot(ps...; layout = length(ps))
end

## ----------------------------------------------------------------------------
function plot_poldist!(p, M::SimModel; 
        hits_count = 3000,
        maxiters = 100 * hits_count,
        static_th = 0.8, rand_th = 0.1,
        pkwargs = Dict(), skwargs = Dict()
    )
    
    vatp_range, vg_ranges = vatpvg_ranges(M)
    vgL, vgU = minimum(first.(vg_ranges)), maximum(last.(vg_ranges))
    Δvg = step(first(vg_ranges))
    vg_range = vgL:Δvg:vgU
    mX, MX = lXgamma(M)
    vatpvgN = sum(length.(vg_ranges))

    hits = 0
    iters = 0
    vatps, vgs, mss = [], [], []

    while vatpvgN > 1 && hits < hits_count && iters < maxiters
        iters += 1
        vatp = rand(vatp_range)
        vg = rand(vg_range)
        !haskey(M.Xb, vatp) && continue
        !haskey(M.Xb[vatp], vg) && continue

        lX = M.Xb[vatp][vg]
        lX/MX < static_th && rand_th > rand() && continue
        
        ms = 10.0 * lX/MX
        push!(vatps, vatp); push!(vgs, vg); push!(mss, ms)
        hits += 1
    end

    scatter!(p, vatps, vgs; 
        color = :black, label = "", ms = mss, 
        alpha = 0.2, skwargs...
    )
    plot!(p; pkwargs...)
end
plot_poldist(M::SimModel; kwargs...) = plot_poldist!(plot(xlabel = "vatp", ylabel = "vg"), M;  kwargs...)


## ----------------------------------------------------------------------------
function plot_polgrid!(p, M; divnum = 50, sparams...)
    def_sparams = (;label = "", color = :blue, alpha = 0.3, lw = 1)
    # vatp
    let 
        vatp_range, vg_ranges = vatpvg_ranges(M)
        sep = length(vatp_range) < divnum ? 1 : length(vatp_range) ÷ divnum
        vatpis = 1:sep:length(vatp_range)
        for vatpi in vatpis
            vatp = vatp_range[vatpi]
            vg_range = vg_ranges[vatpi]
            isempty(vg_range) && continue
            vgL, vgU = minmax(vg_range)
            plot!(p, [vatp, vatp], [vgL, vgU]; def_sparams..., sparams...)
        end
    end

    # vg
    let
        vg_range, vatp_ranges = vgvatp_ranges(M)
        sep = length(vg_range) < divnum ? 1 : length(vg_range) ÷ divnum
        @show sep length(vg_range)
        vgis = 1:sep:length(vg_range)
        for vgi in vgis
            vg = vg_range[vgi]
            vatp_range = vatp_ranges[vgi]
            isempty(vatp_range) && continue
            vatpL, vatpU = minmax(vatp_range)
            plot!(p, [vatpL, vatpU], [vg, vg]; def_sparams..., sparams...)
        end
    end
    return p
end

## ----------------------------------------------------------------------------
function plot_polborder!(p, M; divnum = 50, sparams...)
    def_sparams = (;label = "", color = :gray, alpha = 0.5, lw = 10)
    # vatp
    vatp_range, vg_ranges = vatpvg_ranges(M)
    vatps, vgLs, vgUs = [], [], []
    sep = length(vatp_range) < divnum ? 1 : length(vatp_range) ÷ divnum
    vatpis = 1:sep:length(vatp_range)
    for vatpi in vatpis
        vatp = vatp_range[vatpi]
        vg_range = vg_ranges[vatpi]
        isempty(vg_range) && continue
        vgL, vgU = minmax(vg_range)
        push!(vatps, vatp); push!(vgLs, vgL); push!(vgUs, vgU)
    end
    plot!(p, vatps, vgLs; def_sparams..., sparams...)
    plot!(p, vatps, vgUs; def_sparams..., sparams...)

    return p
end