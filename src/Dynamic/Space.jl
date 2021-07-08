struct Space{T}
    net::MetNet
    box::BoxGrid
    vec::Vector{Vector{T}}
    sortperms::Vector{Vector{Int}}

    function Space{elT}(net::MetNet, frees, δ::Int; 
            stoitolf = 0.01, 
            nthrs = nthreads(),
        ) where {elT  <: AbstractFloat}

        freeis = rxnindex(net, frees)
        box = BoxGrid(net, freeis, δ)
        tol = Δv(net, freeis) .* stoitolf

        boxlen = length(box)

        ch = UtilsJL.GeneralUtils.chunkedChannel(box; 
            chnklen = max(10_000, div(boxlen, 100)), 
            buffsize = nthrs
        )
        
        fea_count_estimate = 0
        
        prog = Progress(boxlen; dt = 0.5, desc = "Collecting Space... ")
        showvalues() = [
            ("fea_count_estimate", fea_count_estimate), 
            ("boxlen", boxlen), 
            ("polV/boxV", fea_count_estimate/boxlen),
            ("thid", threadid()),
        ]

        subspaces_pool = [Container{Vector{elT}}(10_000) for th in 1:nthrs]

        @threads for _ in 1:(nthrs + 1)
            thid = threadid()
            (thid == 1) && continue

            subspace_th = subspaces_pool[thid]
            net_th = deepcopy(net)
            
            for box_chnk in ch
                for f64v in box_chnk
                    if is_feasible(net_th, freeis, f64v; tol)
                        v = collect(elT, f64v)
                        push!(subspace_th, v)
                        fea_count_estimate += 1 
                    end

                    next!(prog; showvalues)
                end
            end

            

        end # for _ in _threads
        finish!(prog)
        
        vec = concat!(subspaces_pool...)
        
        # sortperms
        sortperms = map(eachindex(frees)) do free_ci
            sortperm(vec; by = (v) -> v[free_ci])
        end

        return new{elT}(net, box, vec, sortperms)
    end

    Space(net::MetNet, frees, δ::Int; kwargs...) = Space{Float16}(net, frees, δ; kwargs...)

end

## ------------------------------------------------------
function Vi(V::Space, free::Symbol)
    ci = free_ci(V, free)
    map((v) -> v[ci], V)
end

## ------------------------------------------------------
free_ids(V::Space) = V.box.dims_ids
free_ids(V::Space, free) = V.box.dims_ids[free_ci(V, free)]
free_ci(V::Space, freeid::Symbol) = findfirst(isequal(freeid), free_ids(V))
free_ci(V::Space, free::Int) = findfirst(isequal(free), eachindex(free_ids(V)))
free_ci(V::Space, frees) = free_ci.([V], frees)
Base.vec(V::Space) = V.vec

## ------------------------------------------------------
# Iterator Interface
Base.IteratorSize(::Space) = Base.HasLength()
Base.IteratorEltype(::Space{T}) where {T} = Base.HasEltype()
Base.eltype(V::Space{T}) where {T} = Vector{T}
Base.length(V::Space) = length(V.vec)
Base.size(V::Space) = tuple(length(V))
Base.size(V::Space, dim) = size(V)[dim]
Base.iterate(V::Space) = iterate(V.vec)
Base.iterate(V::Space, state) = iterate(V.vec, state)

## ------------------------------------------------------
# Array interface
Base.eachindex(V::Space) = eachindex(V.vec)
Base.firstindex(V::Space) = firstindex(V.vec)
Base.lastindex(V::Space) = lastindex(V.vec)

Base.getindex(V::Space, i) = getindex(V.vec, i)
Base.setindex!(V::Space, v, i) = setindex!(V.vec, v, i)

## ------------------------------------------------------
function Base.show(io::IO, V::Space{T}) where {T}
    println(io, "Space{", T, "}")
    println(io, "frees: [", join(string.(free_ids(V)), ", "), "]")
    println(io, "len: ", length(V))
    println(io, V.net)
end
Base.show(V::Space) = show(stdout, V)

## ------------------------------------------------------
Base.empty!(V::Space) = (empty!(V.vec); empty!(V.sortperms); V)
