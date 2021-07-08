struct Space{T}
    net::MetNet
    box::BoxGrid
    vec::Vector{Vector{T}}
    iter::Base.Iterators.Zip

    function Space{elT}(net::MetNet, freeids::Vector{Symbol}, δ::Int; 
            nthrs = nthreads(), filter::Function = (net, v) -> true
        ) where {elT  <: AbstractFloat}
        
        box = BoxGrid(net, freeids, δ)
        vec = subspace(net, box, freeids; nthrs, filter)
        any(length.(vec) .!= length(first(vec))) && error("All dimensions must have equal length")
        new(net, box, vec, zip(vec...))
    end

    Space(net::MetNet, freeids::Vector{Symbol}, δ::Int; kwargs...) = Space{Float64}(net, freeids, δ; kwargs...)

end

## ------------------------------------------------------
function Vi(V::Space, free::Symbol)
    ci = free_ci(V, free)
    V.vec[ci]
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
Base.length(V::Space) = length(first(V.vec))
Base.size(V::Space) = tuple(length(V))
Base.size(V::Space, dim) = size(V)[dim]
Base.iterate(V::Space) = iterate(V.iter)
Base.iterate(V::Space, state) = iterate(V.iter, state)

## ------------------------------------------------------
# Array interface
Base.eachindex(V::Space) = eachindex(first(V.vec))
Base.firstindex(V::Space) = firstindex(first(V.vec))
Base.lastindex(V::Space) = lastindex(first(V.vec))

Base.getindex(V::Space, i::Int) = getindex.(V.vec, i)
Base.getindex(V::Space, i) = getindex.([V], i)

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
