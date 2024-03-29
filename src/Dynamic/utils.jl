## ------------------------------------------------------
# utils
red_resolution(v, n::Int) = v[begin:max(div(end, n), 1):end]
normalizeP!(P::AbstractVector; tol = 1e-5) = (Z = sum(P); (abs(Z - 1.0) > tol) && (P .= P ./ Z); P)

function chuncks(v, c::Int)
    vlen = length(v)
    clen = max(div(vlen, c) + 1, 1)
    chnks = 
    idx0 = firstindex(v)
    while true
        idx1 = min(vlen, idx0 + clen)
        push!(chnks, view(v, idx0:idx1))
        (idx1 >= vlen) && return chnks
        idx0 = idx1 + 1
    end
end

function hist(Vi, P)
    hist = Dict{Float64, Float64}()
    for (v, p) in zip(Vi, P)
        get!(hist, v, 0.0)
        hist[v] += p
    end
    hist
end

function n_ave_conv!(v, n; aux = copy(v))
    idx0 = firstindex(v)
    idx1 = lastindex(v)
    for i in eachindex(v)
        ran = max(idx0, i - n):min(idx1, i + n)
        v[i] = sum(view(aux, ran)) / length(ran)
    end
    return v
end

function tail_ave(v, n)
    idx0 = firstindex(v)
    idx1 = lastindex(v)
    tail = view(v, max(idx0, idx1 - n):idx0)
    mean(tail)
end


function check_stst(w, th, vs...)
    for v in vs
        idx1 = lastindex(v)
        (idx1 <= w) && return false
        ran = (idx1 - w):idx1
        vi = view(v, ran)
        rel_var = abs(std(vi) / mean(vi))
        (rel_var > th) && return false
    end
    return true
end

ismul(v, m) = iszero(rem(v, m))

function findapproxi(v, r::AbstractRange)
    v0, v1 = first(r), last(r)
    (v0 > v1) && error("Only increasing ranges are allowed")
    idx = Int(div(v - v0, step(r), RoundNearest)) + 1
    clamp(idx, firstindex(r), lastindex(r))
end

findapproxv(v, r::AbstractRange) = r[findapproxi(v, r)]