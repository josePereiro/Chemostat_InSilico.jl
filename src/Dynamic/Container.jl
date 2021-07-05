## ------------------------------------------------------
mutable struct Container{T}
    v::Vector{T}
    jump_size::Int
    idx::Int

    function Container{T}(jump_size) where {T}
        @assert jump_size > 0
        new{T}(Vector{T}(undef, jump_size), jump_size, 1)
    end
end

function Base.push!(cont::Container{T}, elm::T) where {T}
    vlen = length(cont.v)
    if cont.idx <= vlen
        cont.v[cont.idx] = elm
        cont.idx += 1
    else
        resize!(cont.v, vlen + cont.jump_size)
        push!(cont, elm)
    end
end

concat!(conts::Container...) = vcat(vec!.(conts)...)

vec!(cont::Container) = (resize!(cont.v, cont.idx - 1); cont.v)