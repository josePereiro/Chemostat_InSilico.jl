import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

import Chemostat_InSilico
const InCh = Chemostat_InSilico
const InLP = InCh.LP_Implement
const InU = InCh.Utilities

using ProgressMeter
using Plots
using Plots.PlotMeasures
using InteractiveUtils

# import GR
# GR.inline("png")

import UtilsJL
const UJL = UtilsJL
using Base.Threads
using Serialization
using Statistics
using Dates

## ----------------------------------------------------------------------------
# Meta
fileid = "2.1"
fig_path(fname) = joinpath(InLP.CH3_FIGURES_DIR, fname)

## ----------------------------------------------------------------------------
# DAT [[:DyM, :TS, :M], Vl, D, ϵ]
DAT = UJL.load_data(InCh.CH3_DAT_BUNDLE_FILE)
varinfo(Main, r"DAT")

## ----------------------------------------------------------------------------
# Find beta
let
    Vl, D, ϵ = DAT[:Vls][1], DAT[:Ds][5], DAT[:ϵs][1]
    M = DAT[:M, Vl, D, ϵ]
    DyMs = DAT[:DyM, Vl, D, ϵ]
    
    # Maxent
    δ = 0.05 # discretization factor
    y = abs(M.net.S[2, 6]) # atp/biomass yield
    maxentf(beta) = (vatp, vg) -> exp(beta*vatp/y) + 0.003
    
    # # Gradient descent
    biom_ider = "biom"
    target = InLP.marginal_av(DyMs[biom_ider]) # dynamic mean
    x0 = 1e2
    x1 = 2e2
    maxΔ = 3e2
    th = 1e-3
    gd_xs = []
    gd_ys = []
    maxiters = 15
    beta0 = InU.grad_desc(;target, x0, x1, maxΔ, th, maxiters) do beta
        MEMs = InLP.get_marginals(maxentf(beta), M, [biom_ider]; δ, verbose = false)
        f = InLP.marginal_av(MEMs[biom_ider])
    end
    
    MEMs = InLP.get_marginals(maxentf(beta0), M, [biom_ider]; δ, verbose = false)

    p1 = InLP.plot_marginals(DyMs, [biom_ider])
    p2 = InLP.plot_marginals(MEMs, [biom_ider])
    plot(p1, p2)
end

## ----------------------------------------------------------------------------
let
    f(x) = x^3 - 10.0 * x + 24

    target = 1.0
    x0 = 0.0
    x1 = 0.01
    gd_xs = []
    gd_ys = []
    xtarget = InU.grad_desc(;target, x0, x1, maxΔ = 1e-4, th = 1e-3) do x
        y = f(x)
        push!(gd_xs, x); push!(gd_ys, y)
        y
    end

    p = plot()
    scatter!(p, gd_xs, gd_ys; label = "gd", alpha = 0.7)
    plot!(p, gd_xs, f.(gd_xs); label = "real", ls = :solid, alpha = 0.7)
    p
end
## ----------------------------------------------------------------------------
# _vs_time_vs_ϵ_vs_D
let
    
    marginf = 0.2
    f = identity
    Ds = DAT[:Ds][1:6]
    Vl = DAT[:Vls] |> first

    for field in [:X_ts, :sg_ts, :sl_ts]
        
        ylabel = replace(string(field), "_ts" => "")
        ps = map(Ds) do D0
        
            vals = getfield.(DAT[:TS, Vl, D0, DAT[:ϵs]], field) 
            ylim = InLP.lims(marginf, vals...)
            
            p = plot(;xlabel = "time", ylabel, 
                title = string("D: ", UJL.sci(D0)))
            for (ϵ, val) in zip(DAT[:ϵs], vals) |> collect |> reverse
                plot!(p, f.(val); 
                    label = ϵ, lw = 4, alpha = 0.7, color = Gray(ϵ * 0.8))
            end
            p
        end
        p0 = plot(ps...; layout = length(ps), legend = false, titlefont = 12)
        
        # saving
        pname = UJL.mysavename("$(ylabel)_vs_time_vs_ϵ_vs_D", "png"; Vl)
        fname = fig_path(string(fileid, "_", pname))
        savefig(p0, fname)
        @info "Plotting" fname
    end
end

## ----------------------------------------------------------------------------
# Steady state_vs_D 
let
    f = identity
    Ds = DAT[:Ds][1:8]
    Vl = DAT[:Vls] |> first

    for field in [:X, :sg, :sl]
        ylabel = field
        p = plot(;xlabel = "D", ylabel, title = "Steady State")
        for ϵ in DAT[:ϵs] |> reverse
            Xs = getfield.(DAT[:M, Vl, Ds, ϵ], field) 
            plot!(p, Ds, f.(Xs .+ 1e-8); label = "", lw = 4, alpha = 0.7, color = Gray(ϵ * 0.8))
        end

        # saving
        pname = UJL.mysavename("steady_state_$(field)_vs_D_vs_ϵ", "png"; Vl)
        fname = fig_path(string(fileid, "_", pname))
        savefig(p, fname)
        @info "Plotting" fname
    end
end

## ----------------------------------------------------------------------------
function soft(v; w = 3)
    v = zero(v)
    for i in eachindex(v)

        li = i - w < firstindex(v) ? i : i - w
        ui = i + w > lastindex(v) ? i : i + w

        v[i] = mean(v[li:ui])
    end
    v
end

## ----------------------------------------------------------------------------
# marginals
let
    
    f = indentity
    Ds =  DAT[:Ds][1:6]
    Vl = DAT[:Vls] |> first
    
    write_lock = ReentrantLock()
    vg_plots = Vector(undef, length(Ds))
    vatp_plots = Vector(undef, length(Ds))
    @threads for (Di, D) in Ds |> enumerate |> collect

        vatp_plot = plot(xlabel = "vatp", ylabel = "pdf", title = string("D: ", UJL.sci(D)))
        vg_plot = plot(xlabel = "vg", ylabel = "pdf", title = string("D: ", UJL.sci(D)))

        for ϵ in DAT[:ϵs]
            M = DAT[:M, Vl, D, ϵ]
            vatp_range, vg_ranges = InLP.vatpvg_ranges(M)

            # collect
            vatp_hist = zeros(length(vatp_range))
            vg_hist = Dict{Float64, Float64}()
            for (ivatp, vatp) in enumerate(vatp_range)
                vg_range = vg_ranges[ivatp]
                for vg in vg_range
                    X = M.Xb[vatp][vg]
                    vatp_hist[ivatp] += X
                    get!(vg_hist, vg, 0.0)
                    vg_hist[vg] += X
                end
            end
            vg_range = vg_hist |> keys |> collect |> sort
            vg_hist = [vg_hist[vg] for vg in vg_range]

            if D == 0.0 # is a delta
                vatp_hist = [1.0; zeros(length(vatp_hist))]
                vatp_range = [first(vatp_range); vatp_range]
            end

            # marginals
            lparams = (;label = "", lw = 4, alpha = 0.7, color =  Gray(ϵ * 0.8))
            plot!(vatp_plot, vatp_range, f(vatp_hist ./ M.X); lparams...)
            plot!(vg_plot, vg_range, f(vg_hist ./ M.X); lparams...)
        end

        lock(write_lock) do
            vatp_plots[Di] = vatp_plot
            vg_plots[Di] = vg_plot
            @info "Done" threadid() Di D
        end
        
    end # for D in Ds

    
    # saving
    params = (;legend = false, titlefont = 10, axistitle = 10)
    p = plot(vatp_plots...; layout = length(vatp_plots), params...)
    pname = "vatp_marginals_vs_D_vs_ϵ.png" 
    fname = fig_path(string(fileid, "_", pname))
    savefig(p, fname)
    @info "Plotting" fname

    p = plot(vg_plots...; layout = length(vg_plots), params...)
    pname = "vg_marginals_vs_D_vs_ϵ.png" 
    fname = fig_path(string(fileid, "_", pname))
    savefig(p, fname)
    @info "Plotting" fname
end

## ----------------------------------------------------------------------------
# polytope_ϵ_animation
let
    Ds = DAT[:Ds]
    ϵs = DAT[:ϵs]
    Vl = DAT[:Vls] |> first

    write_lock = ReentrantLock()

    # @threads 
    for D in Ds |> collect
        
        lock(write_lock) do
            @info "Doing" D threadid() now()
        end

        # box
        marginf = 0.1
        vatpL, vatpU = Inf, -Inf
        vgL, vgU = Inf, -Inf
        for ϵ in ϵs
            M = DAT[:M, Vl, D, ϵ]
            vatpL, vatpU, vgL, vgU = InLP.pol_box(M)
            vatp_margin = abs(vatpL - vatpU) * marginf
            vg_margin = abs(vgL - vgU) * marginf
            vatpL = min(vatpL, vatpL - vatp_margin)
            vatpU = max(vatpU, vatpU + vatp_margin)
            vgL = min(vgL, vgL - vg_margin)
            vgU = max(vgU, vgU + vg_margin)
        end

        # plots
        lock(write_lock) do
            @info "Plotting" D threadid() now()
        end
        params = (;legend = false)
        local ps = []
        for ϵ in ϵs
            M = DAT[:M, Vl, D, ϵ]
            
            p = InLP.plot_politope(M) 
            plot!(p; title = string("D: ", UJL.sci(D), " ϵ: ", UJL.sci(ϵ)), params...)
            
            plot!(p; xlim = [vatpL, vatpU], ylim = [vgL, vgU]) 
            push!(ps, p)
            
        end
        rps = [ps; reverse(ps)]

        # gif
        lock(write_lock) do
            @info "Making Gif" D threadid() now()
        end
        pname = UJL.mysavename("polytope_ϵ_animation", "gif"; D)
        fname = fig_path(string(fileid, "_", pname))
        UJL.save_gif(rps, fname; fps = 2.0)
        lock(write_lock) do
            @info "Done!!" D threadid() fname now()
        end
        
    end # for D in Ds
end

