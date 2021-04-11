## ----------------------------------------------------------------------------
# ϵ scale
let
    p = plot(;title = "ϵ color legend", ylabe = "ϵ")
    for ϵ in INDEX[:ϵs]
        plot!(p, fill(ϵ, 10); lw = 8, color = Gray(ϵ * 0.8), label = ϵ)
    end
    mysavefig(p, "ϵ_legend")
end
