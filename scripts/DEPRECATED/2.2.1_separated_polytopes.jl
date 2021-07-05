## ----------------------------------------------------------------------------
# Separated Polytopes
let
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    Ds =  INDEX[:Ds][1:3:end][1:4]
    ϵs = INDEX[:ϵs]

    ps = Plots.Plot[]

    # Find bigger polytope
    G = (;vgU = -Inf)
    for ϵ in ϵs, D in Ds
        status = idxdat([:status], Vl, D, ϵ, τ)
        status != :stst && continue
        M = idxdat([:M], Vl, D, ϵ, τ)
        vatpL, vatpU, vgL, vgU = Dyn.pol_box(M)
        G.vgU < vgU && (G = (;D, vgU, vatpU))
    end

    for ϵ in ϵs, D in Ds
        M = idxdat([:M], Vl, D, ϵ, τ)
        status = idxdat([:status], Vl, D, ϵ, τ)
        status != :stst && continue

        @info("Doing",(Vl, D, ϵ, τ)); println()
        p = plot(;title = "polytope, D = $(UJL.sci(D)), ϵ = $(UJL.sci(ϵ))", 
            xlabel = "vatp", ylabel = "vg"
        )
        skwargs = (;color = :blue, alpha = 0.3, 
            xlim = [-G.vatpU * 0.1, G.vatpU * 1.1],
            ylim = [-G.vgU * 0.1, G.vgU * 1.1], 
        )
        Dyn.plot_polborder!(p, M)
        Dyn.plot_poldist!(p, M; skwargs)

        # Dynamic marginal
        δ = 0.08
        LP_cache = Dyn.vgvatp_cache(M)
        f(vatp, vg) = M.Xb[vatp][vg] / M.X
        DyMs = Dyn.get_marginals(f, M; δ, LP_cache, verbose = false)
        vatp_av = Dyn.av(DyMs["vatp"]) 
        vg_av = Dyn.av(DyMs["gt"]) 

        kwargs = (;alpha = 0.5, label = "", color = :black)
        hline!(p, [vg_av]; lw = 5, ls = :dash, kwargs...)
        vline!(p, [vatp_av]; lw = 5, ls = :dash, kwargs...)
        # scatter!(p, [vatp_av], [vg_av]; ms = 8, kwargs...)

        push!(ps, p)
    end
    layout = length(ϵs), length(Ds)
    mysavefig(ps, "separated_polytopes"; layout)
end