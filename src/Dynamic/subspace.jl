# ------------------------------------------------------
_sep_freeids(freeid, freeids...) = (freeid, freeids)

function _fva_range(net, box, freeid)
    lb, ub = fva(net, freeid)
    ran = dim_range(box, freeid)
    lbi = findapproxi(lb, ran)
    ubi = findapproxi(ub, ran)
    return ran[lbi:ubi]
end

# remove unfeasibles and trivial
_is_feasible(net) = any(fba(net) .!= 0.0)
function _fill_subspace!(containers::Vector, 
        net::MetNet, box, left_freeids, head;
        filter::Function
    )
    
    curr_freeid, left_freeids = _sep_freeids(left_freeids...)
    curr_free_range = _fva_range(net, box, curr_freeid)

    if isempty(left_freeids)
        # base case
        for cur_freev in curr_free_range
            freev = [head; cur_freev]
            !filter(net, freev) && continue
            # all(iszero.(freev)) && continue # remove trivial
            for (vi, containeri) in zip(freev, containers)
                push!(containeri, vi)
            end
        end
    else
        # recursive step
        for cur_freev in curr_free_range
            _head = [head; cur_freev]
            fixxing(net, curr_freeid, cur_freev) do
                _fill_subspace!(containers, net, box, left_freeids, _head; filter)
            end
        end
    end
end

# one container per dim
_build_containers(freeids, jump_size = 10_000) = [Container{Float64}(jump_size) for free in freeids]

function _containers_pool_size(freeids, containerss_pool)
    lens = zeros(Int, length(freeids))
    for containers in containerss_pool
        for i in eachindex(freeids)
            lens[i] += length(containers[i])
        end
    end
    return lens
end

function _reduce_containers_pool!(containerss_pool)

    # one per thread
    acc_containers = Dict()
    for containers in containerss_pool
        # one vector per dimension
        for i in eachindex(containers)
            !haskey(acc_containers, i) && (acc_containers[i] = containers[i]; continue)
            push!(acc_containers[i], vec!(containers[i])...)
            empty!(containers[i])
        end
    end
    map(1:length(acc_containers)) do i
        vec!(acc_containers[i])
    end
end

function subspace(net::MetNet, box::BoxGrid, freeids::Vector{Symbol};
        # stoitof = 0.01, 
        filter::Function = (net, v) -> true,
        nthrs::Int = nthreads()
    )

    # Base
    if length(freeids) == 1
        _fill_subspace!(_build_containers(freeids), net, box, freeids, Float64[]; filter)
    end

    first_freeid, left_freeids = _sep_freeids(freeids...)
    first_free_range = _fva_range(net, box, first_freeid)

    # recursive step
    containerss_pool = [_build_containers(freeids) for _ in 1:nthrs]
    nets_pool = [deepcopy(net) for _ in 1:nthrs]
    
    prog = Progress(length(first_free_range); dt = 0.5, desc = "Collecting... ")
    function showvalues() 
        estimation = sum(_containers_pool_size(freeids, containerss_pool))
        
        return [
            ("polV estimation", estimation), 
            ("thid", threadid())
        ]
    end

    @threads for cur_freev in first_free_range
        
        thid = threadid()
        containers = containerss_pool[thid]
        net_th = nets_pool[thid]
        head = [cur_freev]

        fixxing(net_th, first_freeid, cur_freev) do
            _fill_subspace!(containers, net_th, box, left_freeids, head; filter)
        end
        next!(prog; showvalues)
    end
    finish!(prog)
    
    _reduce_containers_pool!(containerss_pool)
end