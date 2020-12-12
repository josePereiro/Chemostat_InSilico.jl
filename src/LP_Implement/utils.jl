
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
    vg_ranges = Vector{typeof(vatp_range)}(undef, length(vatp_range))
    for (i, vatp) in enumerate(vatp_range)
        vgL, vgU = fixxing(net, vatp_idx, vatp) do 
            fva(net, vg_idx)
        end
        @inbounds vg_ranges[i] = vrange(vgL - vg_margin, vgU + vg_margin, θvg)
    end
    return vatp_range, vg_ranges
end
vatpvg_ranges(M::SimModel; kwargs...) = 
    vatpvg_ranges(M.net, M.θvatp, M.vatp_idx, M.θvg, M.vg_idx; kwargs...)

## ---------------------------------------------------------
function plot_res(M::SimModel; f = (x) -> x)
    p1 = plot(xlabel = "time", ylabel = "conc")
    plot!(p1, f.(M.sg_ts); label = "sg", lw = 3)
    plot!(p1, f.(M.sl_ts); label = "sl", lw = 3)
    
    p2 = plot(xlabel = "time", ylabel = "X")
    plot!(p2, f.(M.X_ts); label = "X", lw = 3)
    
    p = plot([p1, p2]...)
end
    