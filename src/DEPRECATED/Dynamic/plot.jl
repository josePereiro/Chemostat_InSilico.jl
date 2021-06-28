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
    plot_poldist!(p4, M; static_th = 0.8)
    
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
            @error string("Doing $rxn\n", UJL.err_str(err))
        end
    end
    plot(ps...; layout = length(ps))
end

## ----------------------------------------------------------------------------
function plot_poldist!(p, 
        M; static_th = 0.8, ssize = 1000, 
        pkwargs = Dict(), skwargs = Dict(), 
        min_ms = 3.0, max_ms = 15.0
    )

    Xlav = M.X / vatpvgN(M)
    mX, MX = lXgamma(M)
    coords = []
    for (vatp, lXb) in M.Xb
        for (vg, X) in lXb
            X/Xlav < static_th && continue
            push!(coords, (vatp, vg, X/MX))
        end
    end
    shuffle!(coords)

    sample = length(coords) < ssize ? coords : coords[1:ssize]
    xs = getindex.(sample, 1)
    ys = getindex.(sample, 2)
    ms = max.(min_ms, getindex.(sample, 3) .* max_ms)
    scatter!(p, xs, ys;
        color = :black, label = "", 
        ms, alpha = 0.2, skwargs...
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

## ----------------------------------------------------------------------------
function plot_proj(M, ider1 = "biom", ider2 = "gt"; bins = 50, skwargs...)
    p2d = proj2d(M, ider1, ider2; bins)
    p = plot(;xlabel = ider1, ylabel = ider2)
    for (val1, (idx2L, idx2U)) in p2d
        plot!(p, [val1, val1], [idx2L, idx2U]; 
            label = "", color = :black, 
            lw = 8, alpha = 0.3,
            skwargs...
        )
    end
    p
end