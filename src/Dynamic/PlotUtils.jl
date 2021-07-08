# ## ---------------------------------------------------------
# function plot_frees(net::MetNet, frees; 
#         δ = 40,
#         verbose = true, 
#         m = (8, :square),
#         pltkwargs...
#     )
#     freeis = rxnindex(net, frees)
#     frees = rxns(net, freeis)

#     box = BoxGrid(net, freeis, δ)
#     V = length(box)

#     # Start collecting, TODO: make threaded
#     tolf = 0.01
#     tol = Δv(net, freeis) .* tolf

#     boxc = 0
#     polc = 0

#     cid = ("plot_frees", hash(net), δ, hash(freeis))
#     polV = lcache(cid; verbose) do
        
#         polV_ = [Float64[] for id in frees]

#         prog = Progress(V; dt = 0.1, desc = "Collecting... ")
#         for freev in box
#             fullv = fixxing(net, freeis, freev; tol) do
#                 fba(net)
#             end
            
#             # info
#             boxc += 1
#             verbose && next!(prog; showvalues = [(:polfrac, polc/boxc)])

#             # feasible
#             isfea = any(fullv .!= 0.0)
#             !isfea && continue

#             # collect
#             for (v, id) in zip(freev, eachindex(frees))
#                 push!(polV_[id], v)
#             end
            
#             # info
#             polc += 1
            
#         end
#         verbose && finish!(prog)

#         return polV_
#     end

#     # --------------------------------------------
#     println()
#     verbose && @info("Ploting")
#     ps = Plots.Plot[]
#     defaults = (;label = "", m, alpha = 0.5, stroke = false)

#     for i in 1:length(frees), j in (i + 1):length(frees)
#         freex, freey = frees[i], frees[j]
#         polVx, polVy = polV[i], polV[j]
        
#         verbose && @show freex, freey
        
#         p = plot(;xlabel = freex, ylabel = freey)

#         # reduce redundancy
#         pool = Set{Tuple{Float64, Float64}}()
#         curr_vx = nothing
#         for si in sortperm(polVx)
#             vx, vy = polVx[si], polVy[si]

#             if vx != curr_vx # new vx
#                 # plot if not empty
#                 !isempty(pool) && scatter!(p, first.(pool), last.(pool); color = :red,
#                     defaults..., pltkwargs...)

#                 empty!(pool)
#                 curr_vx = vx
#             end
#             push!(pool, (vx, vy))
#         end
        
#         push!(ps, p)
#     end
#     return ps

# end

## ------------------------------------------------------
function marginal_plots(f::Function, V; normalize = true, bins = length(V))
    ps = Plots.Plot[]
    # marginals
    freeids = free_ids(V)
    ridxs = view(shuffle(eachindex(V)), 1:bins)
    for id in freeids
        Vi_ = Vi(V, id)
        Vi_ = view(Vi_, ridxs)
        P = f.(Vi_)
        normalize && normalize!(P)
        p = bar(hist(Vi_, P); label = "", xlabel = string(id), color = :black)
        push!(ps, p)
    end
    ps
end
marginal_plots(V; normalize = true, bins = length(V)) = 
    marginal_plots((v) -> 1.0, V; normalize, bins)

## ------------------------------------------------------
function pol_plots(V; bins = 10000)
    frees = free_ids(V)
    ps = Plots.Plot[]
    ridxs = view(shuffle(eachindex(V)), 1:bins)
    for i1 in 1:length(frees)
        for i2 in (i1 + 1):length(frees)
            id1, id2 = frees[i1], frees[i2]
            V1, V2 = Vi(V, id1), Vi(V, id2)
            p = scatter( view(V1, ridxs), view(V2, ridxs);
                label = "", 
                xlabel = string(id1),
                ylabel = string(id2)
            )
            push!(ps, p)
        end
    end
    ps
end
