## ------------------------------------------------------
mutable struct Container{T}
    v::Vector{T}
    jump_size::Int
    idx::Int

    function Container{T}(jump_size) where {T}
        @assert jump_size > 0
        new{T}(Vector{T}(undef, jump_size), jump_size, 1)
    end
    Container{T}() where {T} = Container{T}(10_000)
end

function Base.push!(cont::Container{T}, elms::T...) where {T}
    for elm in elms
        vlen = length(cont.v)
        if cont.idx <= vlen
            cont.v[cont.idx] = elm
            cont.idx += 1
        else
            resize!(cont.v, vlen + cont.jump_size)
            push!(cont, elm)
        end
    end
    return cont
end

Base.length(cont::Container) = (cont.idx - 1)
Base.firstindex(cont::Container) = firstindex(cont.v)
Base.lastindex(cont::Container) = (cont.idx - 1)
Base.isempty(cont::Container) = (cont.idx == 1)
Base.empty!(cont::Container) = (empty!(cont.v); cont.idx = firstindex(cont); cont)

concat!(conts::Container...) = vcat(vec!.(conts)...)

Base.vec(cont::Container) = view(cont.v, firstindex(cont):lastindex(cont))
Base.resize!(cont::Container) = resize!(cont.v, lastindex(cont))
vec!(cont::Container) = (resize!(cont); cont.v)