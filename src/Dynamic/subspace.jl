# ------------------------------------------------------
_sep_freeis(freei, freeis...) = (freei, freeis)

function _fva_range(net, box, freei)
    lb, ub = fva(net, freei)
    ran = dim_range(box, freei)
    lbi = findapproxi(lb, ran)
    ubi = findapproxi(ub, ran)
    return ran[lbi:ubi]
end

function _fill_subspace!(containers::Vector, 
        net, box, left_freeis, head
    )
    
    curr_freei, left_freeis = _sep_freeis(left_freeis...)
    curr_free_range = _fva_range(net, box, curr_freei)

    if isempty(left_freeis)
        # base case
        for cur_freev in curr_free_range
            freev = [head; cur_freev]
            for (vi, containeri) in zip(freev, containers)
                push!(containeri, vi)
            end
        end
    else
        # recursive step
        for cur_freev in curr_free_range
            _head = [head; cur_freev]
            fixxing(net, curr_freei, cur_freev) do
                _fill_subspace!(containers, net, box, left_freeis, _head)
            end
        end
    end
end

# one container per dim
_build_containers(freeis, jump_size = 10_000) = [Container{Float64}(jump_size) for free in freeis]

function _containers_pool_size(freeis, containerss_pool)
    lens = zeros(Int, length(freeis))
    for containers in containerss_pool
        for i in eachindex(freeis)
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

function subspace(net::MetNet, freeis::Vector{Symbol}, δ::Int; 
        nthrs::Int = nthreads()
    )

    box = BoxGrid(net, freeis, δ)

    # Base
    if length(freeis) == 1
        _fill_subspace!(_build_containers(freeis), net, box, freeis, Float64[])
    end

    first_freei, left_freeis = _sep_freeis(freeis...)
    first_free_range = _fva_range(net, box, first_freei)

    # recursive step
    containerss_pool = [_build_containers(freeis) for _ in 1:nthrs]
    nets_pool = [deepcopy(net) for _ in 1:nthrs]
    
    prog = Progress(length(first_free_range); dt = 0.5, desc = "Collecting... ")
    function showvalues() 
        estimation = sum(_containers_pool_size(freeis, containerss_pool))
        
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

        fixxing(net_th, first_freei, cur_freev) do
            _fill_subspace!(containers, net_th, box, left_freeis, head)
        end
        next!(prog; showvalues)
    end
    finish!(prog)
    
    _reduce_containers_pool!(containerss_pool)
end