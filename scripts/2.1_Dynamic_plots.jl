import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
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
    using Random
end

## ----------------------------------------------------------------------------
# Meta
fileid = "2.1"
fig_path(fname) = joinpath(InLP.DYN_FIGURES_DIR, fname)
Base.first(v::Vector, i) = v[firstindex(v):(min(lastindex(v), i))]
Base.last(v::Vector, i) = v[max(firstindex(v), length(v) - i + 1):lastindex(v)]

## ----------------------------------------------------------------------------
# Load DAT
# INDEX[:DFILE, Vl, D, ϵ, τ]
INDEX = UJL.load_data(InCh.DYN_DATA_INDEX_FILE)

function getdat(dk, indexks...)
    FILE = INDEX[:DFILE, indexks...]
    if FILE isa UJL.ITERABLE
        dat = []
        for F in FILE
            datum = deserialize(F)[dk]
            push!(dat, datum)
        end
        return dat
    else
        dat = deserialize(FILE)
        return dat[dk]
    end
end

## ----------------------------------------------------------------------------
# PLOTS
PS = UJL.DictTree()
function mysavefig(p, pname; params...)
    pname = UJL.mysavename(pname, "png"; params...)
    fname = joinpath(InLP.DYN_FIGURES_DIR, string(fileid, "_", pname))
    PS[pname] = deepcopy(p)
    savefig(p, fname)
    @info "Plotting" fname
end

## ----------------------------------------------------------------------------
# All Polytopes
@time let
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    Ds =  INDEX[:Ds] |> reverse
    colors = Plots.distinguishable_colors(length(Ds))

    for ϵ in INDEX[:ϵs]
        p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg")
        MD = -Inf

        # Find bigger polytope
        G = (;vgub = -Inf)
        for (D, color) in zip(Ds, colors)
            M = getdat(:M, Vl, D, ϵ, τ)
            M.X < INDEX[:death_th] && continue
            vgub = M.net.ub[M.vg_idx]
            G.vgub < vgub && (G = (;D, vgub, color))
        end

        # full polytope
        M = getdat(:M, Vl, G.D, ϵ, τ)
        InLP.plot_politope!(p, M; rand_th = 0.0, 
            skwargs = (;G.color, alpha = 0.4)
        )
        
        # Other polytopes
        for (D, color) in zip(Ds, colors)
            M = getdat(:M, Vl, D, ϵ, τ)
            M.X < INDEX[:death_th] && continue
            @info "Doing" Vl, D, ϵ, τ
            InLP.plot_politope!(p, M; rand_th = 1.0, 
                hits_count = 100, maxiters = 1e4, 
                skwargs = (;color, alpha = 0.9)
            )
            vgub = M.net.ub[M.vg_idx]
            G.vgub < vgub && (G = (;D, vgub, color))
        end
        
        mysavefig(p, "dyn_stst_polytopes"; Vl, ϵ, τ)
    end
end

## ----------------------------------------------------------------------------
# # D vs vatp
# @time let
#     Vl = INDEX[:Vls] |> first
#     τ =  INDEX[:τs] |> first
#     Ds =  INDEX[:Ds]
#     colors = Plots.distinguishable_colors(length(Ds))
#     maxp = 1000
#     static_th = 0.8

#     for ϵ in INDEX[:ϵs]
#         vatpp = plot(;title = "D vs vatp", xlabel = "D", ylabel = "vatp")
#         vgp = plot(;title = "D vs vg", xlabel = "D", ylabel = "vg")
#         MD = -Inf
#         for (D, color) in zip(Ds, colors)
#             M = getdat(:M, Vl, D, ϵ, τ)
#             M.X < INDEX[:death_th] && continue
#             @info "Doing" Vl, D, ϵ, τ

#             vatp_range, vg_ranges = InLP.vatpvg_ranges(M)
#             mX, MX = InLP.lXgamma(M)

#             c = 1
#             for (vatpi, vatp) in vatp_range |> enumerate |> collect |> shuffle!
#                 for vg in vg_ranges[vatpi] |> collect |> shuffle!

#                     c == maxp && break
#                     !haskey(M.Xb, vatp) && continue
#                     !haskey(M.Xb[vatp], vg) && continue

#                     lX = M.Xb[vatp][vg]
#                     if lX/MX > static_th
#                         ms = max(1.0, 10.0 * lX/MX)
#                         scatter!(vatpp, [D], [vatp]; color, label = "", ms, alpha = 0.8)
#                         scatter!(vgp, [D], [vg]; color, label = "", ms, alpha = 0.8)
#                         c+= 1
#                     end
#                 end
#             end
#         end
#         mysavefig(vatpp, "dyn_stst_vatp_vs_D"; Vl, ϵ, τ)
#         mysavefig(vgp, "dyn_stst_vg_vs_D"; Vl, ϵ, τ)
#     end
# end

## ----------------------------------------------------------------------------
# _vs_time_vs_ϵ_vs_D
let
    
    marginf = 0.2
    f = identity
    Ds = INDEX[:Ds][1:3:9]
    fields = [:X_ts, :sg_ts, :sl_ts]
    cparams = (;titlefont = 12, lw = 4, alpha = 0.7)

    for Vl in INDEX[:Vls], τ in INDEX[:τs]
        ps = []
        for field in fields
            for D0 in Ds
                ylabel = replace(string(field), "_ts" => "")
                vals = getfield.(getdat(:TS, Vl, D0, INDEX[:ϵs], τ), field) 
                ylim = InLP.lims(marginf, vals...)
                
                p = plot(;xlabel = "time", ylabel, 
                    title = string("D: ", UJL.sci(D0)))
                for (ϵ, val) in zip(INDEX[:ϵs], vals) |> collect |> reverse
                    plot!(p, f.(val); label = ϵ, color = Gray(ϵ * 0.8), cparams...)
                end
                push!(ps, p)
            end
        end
        
        M, N = length(fields), length(Ds)
        p0 = plot(ps...; layout = grid(M, N), 
            legend = false, size = [400 * M, 300 * N])
        mysavefig(p0, "time_series_vs_ϵ_vs_D"; Vl, τ)
    end
end

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
                Xs = getfield.(getdat(:M, Vl, Ds, ϵ, τ), field) 
                plot!(p, Ds, f.(Xs .+ 1e-8); label = "", lw = 4, alpha = 0.7, color = Gray(ϵ * 0.8))
            end

            # saving
            mysavefig(p, "$(ylabel)_vs_D_vs_ϵ"; Vl, τ)
        end
    end
end

## ----------------------------------------------------------------------------
# marginals
let
    f = identity
    Ds =  INDEX[:Ds][1:3:9]
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    
    write_lock = ReentrantLock()
    vg_plots = Vector(undef, length(Ds))
    vatp_plots = Vector(undef, length(Ds))

    @time for (Di, D) in Ds |> enumerate 

        vatp_plot = plot(xlabel = "vatp", ylabel = "pdf", title = string("D: ", UJL.sci(D)))
        vg_plot = plot(xlabel = "vg", ylabel = "pdf", title = string("D: ", UJL.sci(D)))

        for ϵ in INDEX[:ϵs] |> sort |> reverse
            M = getdat(:M, Vl, D, ϵ, τ)

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
            lparams = (;label = ϵ, lw = 4, alpha = 0.7, color =  Gray(ϵ * 0.8), 
                )
            plot!(vatp_plot, vatp_range, f(vatp_hist ./ M.X); lparams...)
            plot!(vg_plot, vg_range, f(vg_hist ./ M.X); lparams...)
        end

        lock(write_lock) do
            vatp_plots[Di] = vatp_plot
            vg_plots[Di] = vg_plot
            # @info "Done" threadid() Di D
        end
        
    end # for D in Ds
    
    # saving
    params = (;legend = false, titlefont = 10, axistitle = 10)
    M, N = 1, 3
    p = plot(deepcopy.(vatp_plots)...; layout = grid(M, N), 
    size = [400 * N, 300 * M], params...)
    mysavefig(p, "dyn_vatp_marginals_vs_D_vs_ϵ")
    
    p = plot(deepcopy.(vg_plots)...; layout = grid(M, N), 
        size = [400 * N, 300 * M], params...)
    mysavefig(p, "dyn_vg_marginals_vs_D_vs_ϵ")

    ps = [vatp_plots; vg_plots]
    M, N = 2, 3
    p = plot(deepcopy.(ps)...; layout = grid(M, N), 
        size = [400 * N, 300 * M], params...)
    mysavefig(p, "dyn_marginals_vs_D_vs_ϵ")

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
            
#             p = InLP.plot_politope(M) 
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

