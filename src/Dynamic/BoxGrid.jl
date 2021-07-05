
_range(v0, v1, N::Int) = range(v0, v1; length = N)
_range(v0, v1, s::Float64) = range(v0, v1; step = s)
_range(b::Tuple, δ) = _range(b..., δ)

## ------------------------------------------------------
struct BoxGrid
    dims_ids::Vector{Symbol}
    dims_bounds::Vector{Tuple}
    dims_δs::Vector{Union{Int, Float64}}
    dims_ranges::Vector{StepRangeLen}
    dims_steps::Vector{Float64}
    iter::Base.Iterators.ProductIterator
    D::Int

    function BoxGrid(dat::Vector)
        isempty(dat) && error("No dimension info")
        dims_ids = Symbol[]
        dims_bounds = Tuple[]
        dims_δs = Union{Int, Float64}[]
        for (id, lb, ub, δ) in dat
            (lb > ub) && 
                error("Bad bounds, lb >= ub. At $id got ", (lb, ub))
            push!(dims_ids, id)
            push!(dims_bounds, (lb, ub))
            push!(dims_δs, δ)
        end
        dims_ranges = _range.(dims_bounds, dims_δs)
        dims_steps = step.(dims_ranges)
        iter = Iterators.product(dims_ranges...)
        D = length(dims_ids)
        new(
            dims_ids, dims_bounds, dims_δs, 
            dims_ranges, dims_steps, iter, D
        )
    end
end
function Base.show(io::IO, b::BoxGrid) 
    println(io, "BoxGrid(;")
    tab = "    "
    println(io, tab, "# id = (lb, ub, δ)")
    params = String[]
    for (id, (lb, ub), δ) in zip(b.dims_ids, b.dims_bounds, b.dims_δs)
        push!(params, string(tab, id, " = ", (lb, ub, δ)))
    end
    println(io, join(params, ",\n"))
    println(io, ")")
end
Base.show(b::BoxGrid) = show(stdout, b)

## ------------------------------------------------------
# Iterator Interface
Base.IteratorSize(::BoxGrid) = Base.IteratorSize(Base.Iterators.ProductIterator)
Base.IteratorEltype(::BoxGrid) = Base.IteratorEltype(Base.Iterators.ProductIterator)
Base.eltype(b::BoxGrid) = eltype(b.iter)
Base.length(b::BoxGrid) = length(b.iter)
Base.iterate(b::BoxGrid) = iterate(b.iter)
Base.iterate(b::BoxGrid, state) = iterate(b.iter, state)
Base.size(b::BoxGrid) = tuple(length.(b.dims_ranges)...)
Base.size(b::BoxGrid, dim) = size(b)[dim]