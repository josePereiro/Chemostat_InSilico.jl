## ----------------------------------------------------------------------------
let
    # LEGEND
    leg_p = plot(;title = "legend", xaxis = nothing, ylabel = "ϵ")
    for ϵ in MINDEX[:ϵs]
        color = Gray(ϵ * 0.8)
        plot!(leg_p, fill(ϵ, 10); color, label = string("ϵ: ", ϵ), lw = 8)
    end
    mysavefig(leg_p, "eps_legend") 
end