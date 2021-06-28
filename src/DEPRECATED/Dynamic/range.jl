discretize(v, d; mode = RoundDown) = round(v, mode; digits = Int.(d))

## ----------------------------------------------------------------------------
function vrange(L, U, d::Int)
    dx = 10.0^(-d)
    x0 = round(L, RoundUp; digits = d)
    x1 = round(U, RoundDown; digits = d)
    return x0:dx:x1
end

## ----------------------------------------------------------------------------
function vranges(net::MetNet, δ1, v1_idx, v1_margin, δ2, v2_idx, v2_margin)
    v1L, v1U = fva(net, v1_idx)
    v1_range = vrange(v1L - v1_margin, v1U + v1_margin, δ1)
    i_v1_range = collect(enumerate(v1_range))
    v2_ranges = Vector{typeof(v1_range)}(undef, length(v1_range))
    for (v1i, v1) in enumerate(v1_range)
        v2L, v2U = fixxing(net, v1_idx, v1) do 
            fva(net, v2_idx)
        end
        v2_ranges[v1i] = vrange(v2L - v2_margin, v2U + v2_margin, δ2)
    end
    return v1_range, v2_ranges
end

## ----------------------------------------------------------------------------
vatpvg_ranges(net::MetNet, δvatp, vatp_idx, δvg, vg_idx; 
        vatp_margin::Float64 = 0.0,
        vg_margin::Float64 = 0.0,
    ) = vranges(net, δvatp, vatp_idx, vatp_margin, δvg, vg_idx, vg_margin)

vatpvg_ranges(M::SimModel; kwargs...) = 
    vatpvg_ranges(M.net, M.δvatp, M.vatp_idx, M.δvg, M.vg_idx; kwargs...)

vgvatp_ranges(net::MetNet, δvg, vg_idx, δvatp, vatp_idx; 
        vg_margin::Float64 = 0.0,
        vatp_margin::Float64 = 0.0,
    ) = vranges(net, δvg, vg_idx, vg_margin, δvatp, vatp_idx, vatp_margin)

vgvatp_ranges(M::SimModel; kwargs...) = 
    vgvatp_ranges(M.net, M.δvg, M.vg_idx, M.δvatp, M.vatp_idx; kwargs...)


## ----------------------------------------------------------------------------
function ivranges_th(net_pool::Vector{MetNet}, 
        δ1, v1_idx, v1_margin, δ2, v2_idx, v2_margin
    )

    # ranges
    v1L, v1U = fva(first(net_pool), v1_idx)
    v1_range = vrange(v1L - v1_margin, v1U + v1_margin, δ1)
    i_v1_range = collect(enumerate(v1_range))
    v2_ranges = Vector{typeof(v1_range)}(undef, length(v1_range))
    @threads for (v1i, v1) in i_v1_range
        net = net_pool[threadid()]
        v2L, v2U = fixxing(() -> fva(net, v2_idx), net, v1_idx, v1)
        v2_ranges[v1i] = vrange(v2L - v2_margin, v2U + v2_margin, δ2)
    end
    return i_v1_range, v1_range, v2_ranges
end

## ----------------------------------------------------------------------------
ivatpvg_ranges(net_pool::Vector{MetNet}, δvatp, vatp_idx, δvg, vg_idx; 
        vatp_margin::Float64 = 0.0,
        vg_margin::Float64 = 0.0,
    ) = ivranges_th(net_pool, δvatp, vatp_idx, vatp_margin, δvg, vg_idx, vg_margin)


ivatpvg_ranges(M::SimModel, net_pool::Vector{MetNet}; kwargs...) = 
    ivatpvg_ranges(net_pool, M.δvatp, M.vatp_idx, M.δvg, M.vg_idx; kwargs...)


