import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement

    using ProgressMeter
    using Plots
    using Plots.PlotMeasures
    using InteractiveUtils

    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Serialization
    using Statistics
    using Dates
    using Random
    using FileIO
end

## ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
fig_path(fname) = joinpath(InLP.DYN_FIGURES_DIR, fname)
Base.first(v::Vector, i) = v[firstindex(v):(min(lastindex(v), i))]
Base.last(v::Vector, i) = v[max(firstindex(v), length(v) - i + 1):lastindex(v)]

## ----------------------------------------------------------------------------
# Load DAT
# INDEX[:DFILE, Vl, D, ϵ, τ]
INDEX = UJL.load_data(InCh.DYN_DATA_INDEX_FILE)
idxdat(dk, indexks...) = InLP.idxdat(INDEX, dk, indexks...)

## ----------------------------------------------------------------------------
# PLOTS
mysavefig(p, pname; params...) = 
    InLP.mysavefig(p, pname, InLP.DYN_FIGURES_DIR, fileid; params...)

## ----------------------------------------------------------------------------
# Separated Polytopes
let
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    Ds =  INDEX[:Ds][1:3:end][1:4]
    ϵs = INDEX[:ϵs]

    ps = Plots.Plot[]

    # Find bigger polytope
    G = (;vgU = -Inf)
    for ϵ in ϵs, D in Ds
        M = idxdat([:M], Vl, D, ϵ, τ)
        status = idxdat([:status], Vl, D, ϵ, τ)
        status != :stst && continue
        vatpL, vatpU, vgL, vgU = InLP.pol_box(M)
        G.vgU < vgU && (G = (;D, vgU, vatpU))
    end

    for ϵ in ϵs, D in Ds
        M = idxdat([:M], Vl, D, ϵ, τ)
        status = idxdat([:status], Vl, D, ϵ, τ)
        status != :stst && continue

        @info("Doing",(Vl, D, ϵ, τ)); println()
        p = plot(;title = "polytope, D = $(UJL.sci(D)), ϵ = $(UJL.sci(ϵ))", 
            xlabel = "vatp", ylabel = "vg"
        )
        skwargs = (;color = :blue, alpha = 0.3, 
            xlim = [-G.vatpU * 0.1, G.vatpU * 1.1],
            ylim = [-G.vgU * 0.1, G.vgU * 1.1], 
        )
        InLP.plot_polborder!(p, M)
        hits_count = Int(1e3)
        InLP.plot_poldist!(p, M; rand_th = 1.0, static_th = 0.05, 
            hits_count, skwargs, 
        )

        # Dynamic marginal
        δ = 0.08
        LP_cache = InLP.vgvatp_cache(M; marginf = 1.5)
        f(vatp, vg) = M.Xb[vatp][vg] / M.X
        DyMs = InLP.get_marginals(f, M; δ, LP_cache, verbose = false)
        vatp_av = InLP.av(DyMs["vatp"]) 
        vg_av = InLP.av(DyMs["gt"]) 

        kwargs = (;alpha = 0.5, label = "", color = :black)
        hline!(p, [vg_av]; lw = 5, ls = :dash, kwargs...)
        vline!(p, [vatp_av]; lw = 5, ls = :dash, kwargs...)
        # scatter!(p, [vatp_av], [vg_av]; ms = 8, kwargs...)

        push!(ps, p)
    end
    layout = length(ϵs), length(Ds)
    mysavefig(ps, "separated_polytopes"; layout)
    @warn "Exiting"
    exit()
end

## ----------------------------------------------------------------------------
# # Overlaped Polytopes
# @time let
#     Vl = INDEX[:Vls] |> first
#     τ =  INDEX[:τs] |> first
#     Ds =  INDEX[:Ds] |> reverse
#     colors = Plots.distinguishable_colors(length(Ds))

#     for ϵ in INDEX[:ϵs]
#         p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg")
#         MD = -Inf

#         # Find bigger polytope
#         G = (;vgub = -Inf)
#         for (D, color) in zip(Ds, colors)
#             M = idxdat([:M], Vl, D, ϵ, τ)
#             status = idxdat([:status], Vl, D, ϵ, τ)
#             status != :stst && continue
#             vgub = M.net.ub[M.vg_idx]
#             G.vgub < vgub && (G = (;D, vgub, color))
#         end

#         # full polytope
#         M = idxdat([:M], Vl, G.D, ϵ, τ)
#         InLP.plot_poldist!(p, M; rand_th = 0.0, 
#             skwargs = (;G.color, alpha = 0.4)
#         )
        
#         # Other polytopes
#         for (D, color) in zip(Ds, colors)
#             M = idxdat([:M], Vl, D, ϵ, τ)
#             status = idxdat([:status], Vl, D, ϵ, τ)
#             status != :stst && continue
#             @info "Doing" Vl, D, ϵ, τ
#             InLP.plot_poldist!(p, M; rand_th = 1.0, 
#                 hits_count = 100, maxiters = 1e4, 
#                 skwargs = (;color, alpha = 0.9)
#             )
#         end
        
#         mysavefig(p, "dyn_stst_polytopes"; Vl, ϵ, τ)
#     end
# end

# ## ----------------------------------------------------------------------------
# # _vs_time_vs_ϵ_vs_D
# let
    
#     marginf = 0.2
#     f = identity
#     Ds = INDEX[:Ds][1:6:18]
#     fields = [:X_ts, :sg_ts, :sl_ts]
#     cparams = (;lw = 4, alpha = 0.7)

#     for Vl in INDEX[:Vls], τ in INDEX[:τs]
#         ps = Plots.Plot[]
#         for field in fields
#             for D in Ds
#                 ylabel = replace(string(field), "_ts" => "")
#                 vals = getfield.(idxdat([:TS], Vl, D, INDEX[:ϵs], τ), field) 
#                 ylim = InLP.lims(marginf, vals...)
                
#                 p = plot(;xlabel = "time", ylabel, 
#                     title = string("D: ", UJL.sci(D)))
#                 for (ϵ, val) in zip(INDEX[:ϵs], vals) |> collect |> reverse
#                     plot!(p, f.(val); label = ϵ, color = Gray(ϵ * 0.8), cparams...)
#                 end
#                 push!(ps, p)
#             end
#         end
        
#         layout = length(fields), length(Ds)
#         mysavefig(ps, "time_series_vs_ϵ_vs_D"; layout, Vl, τ)
#     end
# end

## ----------------------------------------------------------------------------
# ϵ scale
let
    p = plot(;title = "ϵ color legend", ylabe = "ϵ")
    for ϵ in INDEX[:ϵs]
        plot!(p, fill(ϵ, 10); lw = 8, color = Gray(ϵ * 0.8), label = ϵ)
    end
    mysavefig(p, "ϵ_legend")
end

## ----------------------------------------------------------------------------
# Steady state_vs_D Dynamic
let
    f = identity
    Ds = INDEX[:Ds]

    for Vl in INDEX[:Vls], τ in INDEX[:τs]
        for field in [:X, :sg, :sl]
            ylabel = field
            p = plot(;xlabel = "D", ylabel, title = "Steady State")
            for ϵ in INDEX[:ϵs] |> reverse
                Xs = getfield.(idxdat(:M, Vl, Ds, ϵ, τ), field) 
                plot!(p, Ds, f.(Xs .+ 1e-8); label = "", lw = 4, alpha = 0.7, color = Gray(ϵ * 0.8))
            end

            # saving
            mysavefig(p, "$(ylabel)_vs_D_vs_ϵ"; Vl, τ)
        end
    end
end

## ----------------------------------------------------------------------------
# sg vs eps
let
    Ds =  INDEX[:Ds]
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    ϵs = INDEX[:ϵs]

    p = plot(;title = "dynamic stst", xlabel = "ϵ", ylabel = "sg")
    ps = Plots.Plot[]
    for D in Ds
        xs, ys = [], []
        for ϵ in ϵs |> sort |> reverse
            M = idxdat([:M], Vl, D, ϵ, τ)
            status = idxdat([:status], Vl, D, ϵ, τ)

            @info("Doing", (Vl, D, ϵ, τ), M.X, status)
            status != :stst && continue

            sg = M.sg
            push!(xs, ϵ); push!(ys, sg)
        end
        length(ys) < length(ϵs) - 1 && continue
        plot!(p, xs, ys; label = "", alpha = 0.5, lw = 3)
    end
    mysavefig(p, "eps_vs_sg") 
end

## ----------------------------------------------------------------------------
# marginals
let
    δ = 0.08
    f = identity
    Ds =  INDEX[:Ds]
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    
    write_lock = ReentrantLock()
    vg_plots = Vector{Plots.Plot}(undef, length(Ds))
    vatp_plots = Vector{Plots.Plot}(undef, length(Ds))

    @time for (Di, D) in Ds |> enumerate 

        vatp_plot = plot(xlabel = "vatp", ylabel = "pdf", title = string("D: ", UJL.sci(D)))
        vg_plot = plot(xlabel = "vg", ylabel = "pdf", title = string("D: ", UJL.sci(D)))

        for ϵ in INDEX[:ϵs] |> sort |> reverse
            M = idxdat([:M], Vl, D, ϵ, τ)
            status = idxdat([:status], Vl, D, ϵ, τ)

            @info("Doing", (Vl, D, ϵ, τ), M.X, status)
            status != :stst && continue

            # Dynamic marginal
            LP_cache = InLP.vgvatp_cache(M; marginf = 1.5)
            f(vatp, vg) = M.Xb[vatp][vg] / M.X
            DyMs = InLP.get_marginals(f, M; δ, LP_cache, verbose = false)
            vatp_av = InLP.av(DyMs["vatp"]) 
            vg_av = InLP.av(DyMs["gt"]) 

            # marginals
            lparams = (;label = ϵ, lw = 4, alpha = 0.7, color =  Gray(ϵ * 0.8), legend = false)
            plot!(vatp_plot, DyMs["vatp"]; lparams...)
            vline!(vatp_plot, [vatp_av]; ls = :dash, lparams...)
            plot!(vg_plot, DyMs["gt"]; lparams...)
            vline!(vg_plot, [vg_av]; ls = :dash, lparams...)
        end

        vatp_plots[Di] = vatp_plot
        vg_plots[Di] = vg_plot
        
    end # for D in Ds
    
    # saving
    # layout = 1, length(Ds)
    mysavefig(vatp_plots, "dyn_vatp_marginals_vs_D_vs_ϵ")
    mysavefig(vg_plots, "dyn_vg_marginals_vs_D_vs_ϵ")
    
    ps = Plots.Plot[vatp_plots; vg_plots]
    layout = 2, length(Ds)
    mysavefig(ps, "dyn_marginals_vs_D_vs_ϵ"; layout)

end

## ----------------------------------------------------------------------------
# # polytope_ϵ_animation
# let
#     Ds = DAT[:Ds]
#     ϵs = DAT[:ϵs]
#     Vl = DAT[:Vls] |> first

#     write_lock = ReentrantLock()

#     # @threads 
#     for D in Ds |> collect
        
#         lock(write_lock) do
#             @info "Doing" D threadid() now()
#         end

#         # box
#         marginf = 0.1
#         vatpL, vatpU = Inf, -Inf
#         vgL, vgU = Inf, -Inf
#         for ϵ in ϵs
#             M = DAT[:M, Vl, D, ϵ]
#             vatpL, vatpU, vgL, vgU = InLP.pol_box(M)
#             vatp_margin = abs(vatpL - vatpU) * marginf
#             vg_margin = abs(vgL - vgU) * marginf
#             vatpL = min(vatpL, vatpL - vatp_margin)
#             vatpU = max(vatpU, vatpU + vatp_margin)
#             vgL = min(vgL, vgL - vg_margin)
#             vgU = max(vgU, vgU + vg_margin)
#         end

#         # plots
#         lock(write_lock) do
#             @info "Plotting" D threadid() now()
#         end
#         params = (;legend = false)
#         local ps = []
#         for ϵ in ϵs
#             M = DAT[:M, Vl, D, ϵ]
            
#             p = InLP.plot_poldist(M) 
#             plot!(p; title = string("D: ", UJL.sci(D), " ϵ: ", UJL.sci(ϵ)), params...)
            
#             plot!(p; xlim = [vatpL, vatpU], ylim = [vgL, vgU]) 
#             push!(ps, p)
            
#         end
#         rps = [ps; reverse(ps)]

#         # gif
#         lock(write_lock) do
#             @info "Making Gif" D threadid() now()
#         end
#         pname = UJL.mysavename("polytope_ϵ_animation", "gif"; D)
#         fname = fig_path(string(fileid, "_", pname))
#         UJL.save_gif(rps, fname; fps = 2.0)
#         lock(write_lock) do
#             @info "Done!!" D threadid() fname now()
#         end
        
#     end # for D in Ds
# end

