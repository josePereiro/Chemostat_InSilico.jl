# NOTE: I didn't worked
# - The marginals distribution have artifacts

# ## ------------------------------------------------------
# function _range(lb, ub, ϵ; sigdigits)
#     free_range = range(lb, ub; step = ϵ * 0.7)
#     free_range = map(free_range) do vi
#         # vi = discretize(vi, ϵ)
#         round(vi; sigdigits)
#     end |> unique!
# end

# _sep_freeis(freei, freeis...) = (freei, freeis)

# ## ------------------------------------------------------
# function _fill_subspace!(
#         container::_Container, net, left_freeis, 
#         head_freeis, ϵ, global_tol;
#         sigdigits = 3
#     )
    
#     curr_freei, left_freeis = _sep_freeis(left_freeis...)
#     lb, ub = fva(net, curr_freei)
#     curr_free_range = _range(lb, ub, ϵ[curr_freei]; sigdigits)
#     tol = global_tol[curr_freei]

#     if isempty(left_freeis)
#         # base case
#         head = map(head_freeis) do head_freei
#             lb_, ub_ = fva(net, head_freei)
#             freev = lb_ + ((ub_ - lb_) / 2.0) # average tolerance
#             freev = discretize(freev, ϵ[head_freei])
#             # round(freev; sigdigits)
#         end
#         for cur_freev in curr_free_range
#             _push!(container, [head; cur_freev])
#         end
#     else
#         # recursive step
#         for cur_freev in curr_free_range
#             fixxing(net, curr_freei, cur_freev; tol) do
#                 _fill_subspace!(container, net, left_freeis, head_freeis, ϵ, global_tol; sigdigits)
#             end
#         end
#     end
# end

# function subspace(net, freeis, δ; 
#         stoitolf = 0.0,
#         nthrs = nthreads(),
#         sigdigits = 3
#     )

#     container0 = _Container{Vector{Float64}}(10_000)

#     Δv_ =  Δv(net, freeis)
#     ϵ = Dict(freeis .=> round.(Δv_ ./ δ; sigdigits))
#     stoitol = Δv_ .* stoitolf
#     global_tol = Dict(freeis .=> stoitol)

#     # Base
#     if length(freeis) == 1
#         _fill_subspace!(container0, net, freeis, tuple(), ϵ, global_tol; sigdigits)
#     end

#     first_freei, left_freeis = _sep_freeis(freeis...)
#     lb, ub = fva(net, first_freei)
#     first_free_range = _range(lb, ub, ϵ[first_freei]; sigdigits)
#     tol = global_tol[first_freei]

#     # recursive step
#     head_freeis = freeis[1:end - 1]
#     containers_pool = [container0; [deepcopy(container0) for _ in 2:nthrs]]
#     nets_pool = [deepcopy(net) for _ in 1:nthrs]
    
#     prog = Progress(length(first_free_range); dt = 0.5, desc = "Collecting... ")
#     showvalues() = [
#         ("polV estimation", sum(cont.idx for cont in containers_pool)), 
#         ("thid", threadid())
#     ]

#     @threads for cur_freev in first_free_range
#         thid = threadid()
#         container = containers_pool[thid]
#         net_th = nets_pool[thid]
#         fixxing(net_th, first_freei, cur_freev; tol) do
#             _fill_subspace!(container, net_th, left_freeis, head_freeis, ϵ, global_tol; sigdigits)
#         end
#         next!(prog; showvalues)
#     end
#     finish!(prog)
#     _vcat!(containers_pool)
# end
