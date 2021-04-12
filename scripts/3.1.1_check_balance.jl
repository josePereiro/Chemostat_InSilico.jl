## ----------------------------------------------------------------------------
# Check dyn quality
let

    biom_prods, bioms_drains = [], []
    glc_ins, glc_ups, glc_drains = [], [], []
    # LP_cache = nothing
    for (Vl, D, ϵ, τ) in collect(EXP_PARAMS)
        
        # LOAD
        status = MINDEX[:STATUS, Vl, D, ϵ, τ]
        @info("Doing", (Vl, D, ϵ, τ), status)
        status != :stst && continue
        println()

        DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        M = idxdat([:M0], Vl, D, ϵ, τ)
        
        # COMPUTE
        f(vatp, vg) = M.Xb[vatp][vg] / M.X
        av_z = Dyn.av(DyMs["biom"])
        μ = av_z - M.σ
        growth_bal = (μ - D)/D
        biom_prod = μ * M.X
        bioms_drain = D * M.X

        av_vg = Dyn.av(DyMs["gt"])
        cD_X = M.cg * M.D / M.X
        glc_bal = (-av_vg * M.X + (M.cg - M.sg) * M.D) / (M.cg * M.D)
        glc_in = M.cg * M.D
        glc_drain = M.sg * M.D
        glc_up = av_vg * M.X
        
        # PUSH
        @info("Growth", av_z, μ, M.σ, D, growth_bal)
        @info("GLC uptake", av_vg, cD_X, M.sg, glc_bal)
        println()

        push!(biom_prods, biom_prod)
        push!(bioms_drains, bioms_drain)
        push!(glc_ins, glc_in)
        push!(glc_drains, glc_drain)
        push!(glc_ups, glc_up)
    end

    # PLOT
    biom_bal_p = plot(;title = "Biomass balance", 
        xlabel = "simulation", ylabel = "rate"
    )
    plot!(biom_bal_p, biom_prods; label = "prods", lw = 3)
    plot!(biom_bal_p, bioms_drains; label = "drains", lw = 3)

    biom_corr_p = plot(title = "Biomass balance correlation", 
       xlabel = "production", ylabel = "drain"
    )
    scatter!(biom_corr_p, biom_prods, bioms_drains; label = "", m = 8)
    vals = sort([biom_prods; bioms_drains])
    plot!(biom_corr_p, vals, vals; label = "", color = :black, ls = :dash, alpha = 0.7)
    
    glc_bal_p = plot(;title = "Glc balance", 
        xlabel = "simulation", ylabel = "rate"
    )
    plot!(glc_bal_p, glc_ins; label = "input", lw = 3)
    plot!(glc_bal_p, glc_ups; label = "uptake", lw = 3)
    plot!(glc_bal_p, glc_drains; label = "drain", lw = 3)
    plot!(glc_bal_p, glc_drains .+ glc_ups; label = "up + drain", lw = 3)

    glc_corr_b = plot(title = "Glc balance correlation", 
        xlabel = "production", ylabel = "up + drain"
    )
    scatter!(glc_corr_b,  glc_ins, glc_drains .+ glc_ups; label = "", m = 8)
    vals = sort([glc_ins; glc_drains .+ glc_ups])
    plot!(glc_corr_b, vals, vals; label = "", color = :black, ls = :dash, alpha = 0.7)

    ps = Plots.Plot[biom_bal_p, biom_corr_p, glc_bal_p, glc_corr_b]
    mysavefig(ps, "dyn_balances")

end
