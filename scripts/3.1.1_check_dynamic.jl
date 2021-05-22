## ----------------------------------------------------------------------------
# Check dyn quality
let

    biom_prods, bioms_drains = [], []
    glc_ins, glc_ups, glc_drains = [], [], []
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

# TODO: check marginals for stoi err inconsistency
## ----------------------------------------------------------------------------
# M = idxdat([MEmode, :M], Vl, D, ϵ, τ)
# beta_biom = idxdat([MEmode, :beta_biom], Vl, D, ϵ, τ)
# beta_vg = idxdat([MEmode, :beta_vg], Vl, D, ϵ, τ)
# MEMs = Dyn.get_marginals(M; δ, LP_cache, verbose = false) do vatp_, vg_
#     exp((beta_biom * z(vatp_, vg_)) + (beta_vg * vg(vatp_, vg_)))
# end
# ## ----------------------------------------------------------------------------
# # stoi error
# let
    
#     LP_cache = nothing
#     MEmode = ME_FULL_POLYTOPE
#     for (Vl, D, ϵ, τ) in collect(EXP_PARAMS)
        
#         # LOAD
#         status = MINDEX[:STATUS, Vl, D, ϵ, τ]
#         status != :stst && continue
#         println()

#         M = idxdat([:M0], Vl, D, ϵ, τ)

#         isnothing(LP_cache) && (LP_cache = Dyn.vgvatp_cache(M))

#         δ = 0.08 # marginal discretization factor
#         fX(vatp, vg) = M.Xb[vatp][vg] / M.X
#         DyMs = Dyn.get_marginals(fX, M; δ, LP_cache, verbose = false)
#         # DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        
#         Mavs = [Dyn.av(DyMs[rxn]) for rxn in M.net.rxns]
#         Mavs = [isnan(v) ? 0.0 : v for v in Mavs]
        
#         av_vatp = Dyn.av(DyMs["vatp"])
#         av_vg = Dyn.av(DyMs["gt"])

#         dav_vatp = Dyn.discretize(av_vatp, M.δvatp)
#         dav_vg = Dyn.discretize(av_vg, M.δvg)
#         # @show dav_vatp dav_vg
#         Cavs = LP_cache[dav_vatp][dav_vg]
        
#         Mmaxerr = maximum(abs.(M.net.S * Mavs))
#         Cmaxerr = maximum(abs.(M.net.S * Cavs))
#         @info(
#             "Doing", (Vl, D, ϵ, τ), 
#             status, 
#             av_vatp, av_vg, 
#             Mmaxerr, 
#             s = "------",
#             dav_vatp, dav_vg,
#             Cmaxerr, 
#         ); 

#         sep = 25
#         println("\n", "-"^45)
#         @info("Marginals vs Cache")
#         for (rxn, Mav, Cav) in zip(M.net.rxns, Mavs, Cavs)
#             println(rxn, 
#                 rpad(string(" \tMav: ", Mav), sep), 
#                 rpad(string(" \tCav: ", Cav), sep)
#             )
#         end

#         println("\n", "-"^45)
#         @info("FBA")
#         Dyn.fixxing(M.net, "gt", av_vg) do
#             Dyn.fixxing(M.net, "vatp", av_vatp) do
#                 fbaosol = Dyn.fba(M.net)
#                 for (rxn, v, lb, ub) in zip(M.net.rxns, fbaosol, M.net.lb, M.net.ub)
#                     println(rxn, 
#                         rpad(string(" \tv: ", v), sep),
#                         rpad(string(" \tlb: ", lb), sep),
#                         rpad(string(" \tub: ", ub), sep)
#                     )
#                 end        
#             end
#         end
        
#         println("\n"^3)
#     end
# end

# ## ----------------------------------------------------------------------------
# let
#     _discretize(n, Δ) = round(n / Δ) * Δ
#     _discretize(1.0, 3)
# end
# ## ----------------------------------------------------------------------------
# # stoi error
# let
    
#     LP_cache = nothing
#     MEmode = ME_FULL_POLYTOPE
#     for (Vl, D, ϵ, τ) in collect(EXP_PARAMS)
        
#         # LOAD
#         status = MINDEX[:STATUS, Vl, D, ϵ, τ]
#         status != :stst && continue
#         println()

        
#         M = idxdat([:M0], Vl, D, ϵ, τ)
#         isnothing(LP_cache) && (LP_cache = Dyn.vgvatp_cache(M))

#         δ = 0.08 # marginal discretization factor
#         fX(vatp, vg) = M.Xb[vatp][vg] / M.X
#         DyMs = Dyn.get_marginals(fX, M; δ, LP_cache, verbose = true)
#         # DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        
#         Mavs = [Dyn.av(DyMs[rxn]) for rxn in M.net.rxns]
#         Mavs = [isnan(v) ? 0.0 : v for v in Mavs]

#         av_vatp = Dyn.av(DyMs["vatp"])
#         av_vg = Dyn.av(DyMs["gt"])

#         dav_vatp = Dyn.discretize(av_vatp, M.δvatp)
#         dav_vg = Dyn.discretize(av_vg, M.δvg)
#         Cavs = LP_cache[dav_vatp][dav_vg]

        
#         Mmaxerr = maximum(abs.(M.net.S * Mavs))
#         Cmaxerr = maximum(abs.(M.net.S * Cavs))
        
#         @info(
#             "Doing", (Vl, D, ϵ, τ), 
#             status, 
#             av_vatp, av_vg, 
#             Mmaxerr, 
#             s = "------",
#             dav_vatp, dav_vg,
#             Cmaxerr, 
#         ); println()
#         return
#     end
# end