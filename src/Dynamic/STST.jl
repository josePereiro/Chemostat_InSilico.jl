# ---------------------------------------------------------------
struct STST
    buffs::Vector{Vector{Float64}}
    nbuffs::Int
    bufflen::Int
    absmax::Vector
    state::Vector
    function STST(nbuffs, bufflen)
        @assert nbuffs > 0
        @assert bufflen > 0
        buffs = [zeros(bufflen) for _ in 1:nbuffs]
        absmax = zeros(bufflen)
        curridx = 1
        firstfill = 0 # false
        state = [curridx, firstfill]
        new(buffs, nbuffs, bufflen, absmax, state)
    end
end

function Base.push!(stst::STST, vs...)
    curridx, firstfill = stst.state

    for i in eachindex(stst.buffs)
        stst.buffs[i][curridx] = vs[i]
        stst.absmax[i] = max(abs(vs[i]), stst.absmax[i])
    end

    curridx = (curridx >= stst.bufflen) ? 
        (firstfill = 1) : (curridx + 1)
    stst.state .= [curridx, firstfill]
end

function is_steady(stst, ststth)
    _ , firstfill = stst.state
    (firstfill == 0) && return false

    for i in eachindex(stst.buffs)
        max_ = maximum(stst.buffs[i])
        min_ = minimum(stst.buffs[i])
        norm_ = stst.absmax[i]
        norm_ = iszero(norm_) ? 1.0 : norm_
        var_ = (max_ - min_) / norm_
        (var_ > ststth) && return false
    end
    return true
end
