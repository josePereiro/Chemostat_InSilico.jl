## ------------------------------------------------------
struct PMF{D, P}
    domain::Vector{D}
    prob::Vector{P}
    iter::Base.Iterators.Zip{Tuple{Vector{D}, Vector{P}}}
end
PMF(domain::Vector{D}, prob::Vector{P}) where {D, P} =
    PMF{D, P}(domain, prob, zip(domain, prob))

## ------------------------------------------------------
function _pmf(f::Function, V::Space, free::Symbol; normalize = true)
    ci = free_ci(V, free)

    first_vi = first(V)[ci]
    estimated_len = maximum(size(V.box))
    domain = Container{typeof(first_vi)}(estimated_len)
    prob = Container{Float64}(estimated_len)

    last_vi = first_vi
    acc_prob = 0.0
    for i in V.sortperms[ci]
        v = V[i]
        vi = v[ci]
        
        if (last_vi != vi)
            push!(domain, vi)
            push!(prob, acc_prob)
            acc_prob = 0.0
            last_vi = vi
        end
        acc_prob += f(v)
    end

    domain = vec!(domain)
    prob = vec!(prob)
    
    normalize && (prob ./= sum(prob))

    return PMF(domain, prob)
end

function _pmf(f::Function, V::Space, free; normalize = true)
    ci = free_ci(V, free)
    
    first_vi = first(V)[ci]
    dict = Dict{typeof(first_vi), Float64}()

    for v in V
        vi = v[ci]
        get!(dict, vi, 0.0)
        dict[vi] += f(v)
    end

    domain = collect(keys(dict))
    prob = collect(values(dict))

    normalize && (prob ./= sum(prob))
    return PMF(domain, prob)

end

function _pmf(f::Function, V::Space; normalize = true)
    domain = vec(V)
    prob = map(f, domain)
    normalize && (prob ./= sum(prob))
    return PMF(domain, prob)
end

_is_cannonical_span(V, free) = all(free_ci.([V], free) .== free_ci.([V], free_ids(V)))
function pmf(f::Function, V::Space, free = free_ids(V); normalize = true)
    if _is_cannonical_span(V, free)
        _pmf(f, V; normalize)
    else
        _pmf(f, V, free; normalize)
    end
end

pmf(V::Space, free = free_ids(V); normalize = true) = pmf((v) -> 1.0, V, free; normalize)

## ------------------------------------------------------
# Iterator Interface
Base.IteratorSize(P::PMF) = HasLength(P.iter)
Base.IteratorEltype(P::PMF{T}) where {T} = Base.IteratorEltype(P.iter)
Base.eltype(P::PMF{T}) where {T} = eltype(P.ter)
Base.length(P::PMF) = length(P.iter)
Base.size(P::PMF) = size(P.iter)
Base.size(P::PMF, dim) = size(P)[dim]

Base.iterate(P::PMF) = iterate(P.iter)
Base.iterate(P::PMF, state) = iterate(P.iter, state)

## ------------------------------------------------------
# Array interface
Base.getindex(P::PMF, i) = (getindex(P.domain, i), getindex(P.prob, i))
Base.setindex!(P::PMF, p, i) = setindex!(P.prob, p, i)

Base.eachindex(P::PMF) = eachindex(P.domain)
Base.firstindex(P::PMF) = firstindex(P.domain)
Base.lastindex(P::PMF) = lastindex(P.domain)

## ------------------------------------------------------
# show
function Base.show(io::IO, P::PMF{DT, PT}) where {DT, PT}
    println(io, "PMF{", DT, ", ", PT, "}")
    println(io, "size: ", size(P))
end
Base.show(P::PMF) = show(stdout, P)

## ------------------------------------------------------
function sample(pmf::PMF, d::Int)
    len = length(pmf)
    step = max(div(len, d), 1)
    ran = 1:step:len
    PMF(pmf.domain[ran], pmf.prob[ran])
end

## ------------------------------------------------------
function normalize!(P::PMF; tol = 1e-8)
    Z = sum(P.prob)
    (abs(Z - 1.0) > tol) && (P.prob .= P.prob ./ Z)
    P
end

function _ave(f::Function, P::PMF, ave_::Vector)
    for i in eachindex(P)
        v, p = P[i]
        ave_ .= ave_ .+ (f(v) .* p) 
    end
    ave_
end

function _ave(f::Function, P::PMF, ave_::AbstractFloat)
    for i in eachindex(P)
        v, p = P[i]
        ave_ = ave_ + (f(v) * p) 
    end
    ave_
end

ave(f::Function, P::PMF) = _ave(f, P, zero(f(first(P.domain))))
ave(P::PMF) = ave(identity, P)

var(f::Function, P::PMF; av = ave(f, P)) = 
    ave((v) -> (av .- f(v)).^2, P)
var(P; av = ave(f, P)) = var(identity, P; av)
