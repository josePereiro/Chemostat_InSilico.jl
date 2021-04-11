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
#         Dyn.plot_poldist!(p, M; rand_th = 0.0, 
#             skwargs = (;G.color, alpha = 0.4)
#         )
        
#         # Other polytopes
#         for (D, color) in zip(Ds, colors)
#             M = idxdat([:M], Vl, D, ϵ, τ)
#             status = idxdat([:status], Vl, D, ϵ, τ)
#             status != :stst && continue
#             @info "Doing" Vl, D, ϵ, τ
#             Dyn.plot_poldist!(p, M; rand_th = 1.0, 
#                 hits_count = 100, maxiters = 1e4, 
#                 skwargs = (;color, alpha = 0.9)
#             )
#         end
        
#         mysavefig(p, "dyn_stst_polytopes"; Vl, ϵ, τ)
#     end
# end