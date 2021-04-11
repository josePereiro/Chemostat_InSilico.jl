## ----------------------------------------------------------------------------
# marginals
let
    δ = 0.08
    f = identity
    Ds =  INDEX[:Ds][1:3:end]
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    
    vg_plots = Vector{Plots.Plot}(undef, length(Ds))
    vatp_plots = Vector{Plots.Plot}(undef, length(Ds))

    @time for (Di, D) in Ds |> enumerate 

        vatp_plot = plot(;title = string("D: ", UJL.sci(D)), 
            xlabel = "vatp", ylabel = "pdf"
        )
        vg_plot = plot(;title = string("D: ", UJL.sci(D)), 
            xlabel = "vg", ylabel = "pdf"
        )

        for ϵ in INDEX[:ϵs] |> sort |> reverse
            M = idxdat([:M], Vl, D, ϵ, τ)
            status = idxdat([:status], Vl, D, ϵ, τ)

            @info("Doing", (Vl, D, ϵ, τ), M.X, status)
            status != :stst && continue

            # Dynamic marginal
            LP_cache = Dyn.vgvatp_cache(M)
            f(vatp, vg) = M.Xb[vatp][vg] / M.X
            DyMs = Dyn.get_marginals(f, M; δ, LP_cache, verbose = false)
            vatp_av = Dyn.av(DyMs["vatp"]) 
            vg_av = Dyn.av(DyMs["gt"]) 

            # marginals
            lparams = (;label = ϵ, lw = 4, alpha = 0.7, 
                color =  Gray(ϵ * 0.8), legend = false
            )
            plot!(vatp_plot, DyMs["vatp"]; lparams...)
            vline!(vatp_plot, [vatp_av]; ls = :dash, lparams...)
            plot!(vg_plot, DyMs["gt"]; lparams...)
            vline!(vg_plot, [vg_av]; ls = :dash, lparams...)
        end

        vatp_plots[Di] = vatp_plot
        vg_plots[Di] = vg_plot
        
    end # for D in Ds
    
    # saving
    # layout = 1, length(Ds)
    mysavefig(vatp_plots, "dyn_vatp_marginals_vs_D_vs_ϵ")
    mysavefig(vg_plots, "dyn_vg_marginals_vs_D_vs_ϵ")
    
    ps = Plots.Plot[vatp_plots; vg_plots]
    layout = 2, length(Ds)
    mysavefig(ps, "dyn_marginals_vs_D_vs_ϵ"; layout)

end
