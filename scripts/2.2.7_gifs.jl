# ## ----------------------------------------------------------------------------
# # Progress gifs
# let
#     PROG_FIG_DIR = joinpath(InCh.Dyn.plotsdir(), "progress")
#     UJL.make_group_gif("tslen", PROG_FIG_DIR)
# end

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
#             vatpL, vatpU, vgL, vgU = Dyn.pol_box(M)
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
            
#             p = Dyn.plot_poldist(M) 
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
#         fname = Dyn.plotsdir(string(fileid, "_", pname))
#         UJL.save_gif(rps, fname; fps = 2.0)
#         lock(write_lock) do
#             @info "Done!!" D threadid() fname now()
#         end
        
#     end # for D in Ds
# end

