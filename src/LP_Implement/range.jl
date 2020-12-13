function vrange(L, U, d::Int)
    dx = 10.0^(-d)
    x0 = round(L, RoundUp; digits = d)
    x1 = round(U, RoundDown; digits = d)
    return x0:dx:x1
end

function vatpvg_ranges(net::MetNet, θvatp, vatp_idx, θvg, vg_idx; 
        vatp_margin::Float64 = 0.0,
        vg_margin::Float64 = 0.0,
    )
    # ranges
    vatpL, vatpU = fva(net, vatp_idx)
    vatp_range = vrange(vatpL - vatp_margin, vatpU + vg_margin, θvatp)
    i_vatp_range = collect(enumerate(vatp_range))
    vg_ranges = Vector{typeof(vatp_range)}(undef, length(vatp_range))
    for (vatpi, vatp) in enumerate(vatp_range)
        vgL, vgU = fixxing(net, vatp_idx, vatp) do 
            fva(net, vg_idx)
        end
        vg_ranges[vatpi] = vrange(vgL - vg_margin, vgU + vg_margin, θvg)
    end
    return vatp_range, vg_ranges
end

vatpvg_ranges(M::SimModel; kwargs...) = 
    vatpvg_ranges(M.net, M.θvatp, M.vatp_idx, M.θvg, M.vg_idx; kwargs...)

# TODO finish this threading version
function ivatpvg_ranges(net_pool::Vector{MetNet}, θvatp, vatp_idx, θvg, vg_idx; 
        vatp_margin::Float64 = 0.0,
        vg_margin::Float64 = 0.0,
    )

    # ranges
    vatpL, vatpU = fva(first(net_pool), vatp_idx)
    vatp_range = vrange(vatpL - vatp_margin, vatpU + vg_margin, θvatp)
    i_vatp_range = collect(enumerate(vatp_range))
    vg_ranges = Vector{typeof(vatp_range)}(undef, length(vatp_range))
    @threads for (vatpi, vatp) in i_vatp_range
        tid = threadid()
        net = net_pool[tid]
        vgL, vgU = fixxing(net, vatp_idx, vatp) do 
            fva(net, vg_idx)
        end
        vg_ranges[vatpi] = vrange(vgL - vg_margin, vgU + vg_margin, θvg)
    end
    return i_vatp_range, vatp_range, vg_ranges
end

ivatpvg_ranges(M::SimModel, net_pool::Vector{MetNet}; kwargs...) = 
    ivatpvg_ranges(net_pool, M.θvatp, M.vatp_idx, M.θvg, M.vg_idx; kwargs...)
