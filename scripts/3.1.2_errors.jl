## ----------------------------------------------------------------------------
# Error
let
    Vl = MINDEX[:Vls] |> first
    τ = MINDEX[:τs] |> first
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort

    ps = Plots.Plot[]
    M = -Inf
    for MODsym in ALL_MODELS
        p = plot(;tile = string("Error", MODsym), 
            xlabel = "log ϵ", ylabel = "maximum err"
        )
        xs, ys, yerrs = [], [], []
        for ϵ in ϵs
            errs = []
            @info("Doing", ϵ, MODsym); println()
            for D in Ds
                status = MINDEX[:STATUS, Vl, D, ϵ, τ]
                status != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                
                for rxn in Dyn.RXNS
                    DYN_flx = Dyn.av(DyMs[rxn])
                    ME_flx = Dyn.av(Ms[rxn])
                    (isnan(DYN_flx) || isnan(ME_flx)) && continue

                    err = ((DYN_flx - ME_flx)^2)/abs(DYN_flx)
                    push!(errs, err)
                end
            end

            push!(xs, ϵ)
            max_err = maximum(errs)
            push!(ys, max_err)
            M = max(M, max_err)

        end # for ϵ in ϵs

        noise = xs .* 0.1 .* rand.()
        params = (;alpha = 0.8, color = MOD_COLORS[MODsym])
        plot!(p, log10.(xs .+ noise), ys; 
            label = "", lw = 3, ls = :dash, params...
        )
        scatter!(p, log10.(xs .+ noise), ys; # yerr = yerrs, 
            ms = 8, params..., label = string(MODsym), legend = :topleft
        )
        push!(ps, p)
    end #  for MODsym 

    mysavefig(ps, "eps_vs_err")
    
    for p in ps
        plot!(p; ylim = [0.0, M * 1.1])
    end
    mysavefig(ps, "eps_vs_err_equal_scale")

end