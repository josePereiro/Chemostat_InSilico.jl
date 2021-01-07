function similar_board(Xb, v0)
    T = typeof(v0)
    sXb = Dict{Float64, Dict{Float64, T}}()
    for (vatp, lX) in Xb
        slX = get!(sXb, vatp, Dict{Float64, T}())
        for (vg, _) in lX
            slX[vg] = v0
        end
    end
    sXb
end

function complete_board!(Xb, i_vatp_range, vg_ranges, deflt)
    T = eltype(values(Xb))
    for (vatpi, vatp) in i_vatp_range
        vg_range = vg_ranges[vatpi]
        lXb = haskey(Xb, vatp) ? Xb[vatp] : (Xb[vatp] = T())
        for vg in vg_range
            !haskey(lXb, vg) && (lXb[vg] = deflt)
        end
    end
    Xb
end

function fill_board!(Xb, v0)
    for (vatp, lX) in Xb
        for (vg, _) in lX
            lX[vg] = v0
        end
    end
end