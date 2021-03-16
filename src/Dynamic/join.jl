function get_join!(f, M, join::Dict{Float64, Dict{Float64, Float64}})

    vatp_range, vg_ranges = vatpvg_ranges(M)

    # COLLECTING
    Z = 0.0
    for (vatpi, vatp) in enumerate(vatp_range)
        join_vatp = get!(join, vatp, Dict{Float64, Float64}())
        for vg in vg_ranges[vatpi]
            # call function
            fv = f(vatp, vg)
            join_vatp[vg] = fv
            Z += fv
        end
    end
    @assert Z > 0.0
    @assert Z < Inf

    # NORMALIZING
    for (vatpi, vatp) in enumerate(vatp_range)
        join_vatp = join[vatp]
        for vg in vg_ranges[vatpi]
            join_vatp[vg] /= Z
        end
    end

    return join
end
get_join(f, M) = get_join!(f, M, Dict{Float64, Dict{Float64, Float64}}())
get_join(M) = get_join((args...) -> 1.0, M)

## ----------------------------------------------------------------------------
function ave_over(f::Function, P)
    ave = 0.0
    for (vatp, mP) in P
        for (vg, p) in mP
            ave += f(vatp, vg) * p
        end
    end
    ave
end
