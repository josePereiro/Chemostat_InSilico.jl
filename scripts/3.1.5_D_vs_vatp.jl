## ----------------------------------------------------------------------------
# D vs vatp
let
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort
    τ = MINDEX[:τs] |> first
    Vl = MINDEX[:Vls] |> first

    p = plot(;title = "dynamic stst", xlabel = "D", ylabel = "vatp")
    for ϵ in ϵs
        @info("Doing", ϵ)
        xs, ys = [], []
        for D in Ds
            MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
            DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
            
            vatp = Dyn.av(DyMs["vatp"])
            push!(xs, D) 
            push!(ys, vatp) 
        end
        plot!(p, xs, ys; label = ϵ, alpha = 0.5, lw = 3)
    end
    mysavefig(p, "mu_vs_eps") 
end

## ----------------------------------------------------------------------------
# \BETA vs D vs \EPS
let
    MODsym = ME_FULL_POLYTOPE

    dat_pool = Dict()
    for (Vl, D, ϵ, τ) in EXP_PARAMS
        dat = get!(dat_pool, ϵ, Dict())
        get!(dat, :Ds, [])
        get!(dat, :betas, [])

        status = MINDEX[:STATUS, Vl, D, ϵ, τ] 
        
        @info("Doing", (Vl, D, ϵ, τ), MODsym, status)
        status != :stst && continue
        beta = idxdat([MODsym, :beta_biom], Vl, D, ϵ, τ)
        # beta = rand() # Test

        push!(dat[:Ds], D)
        push!(dat[:betas], beta)
    end

    # PLOTTING D vs BETA
    p = plot(;xlabel = "D", ylabel = "beta")
    for (eps, dat) in dat_pool
        color = Gray(eps * 0.8)
        scatter!(p, dat[:Ds], dat[:betas]; 
            label = "", color, m = 8, alpha = 0.7
        )
    end
    mysavefig(p, "beta_vs_D"; MODsym) 

    # PLOTTING EPS vs BETA
    p = plot(;xlabel = "eps", ylabel = "beta")
    for (eps, dat) in dat_pool
        epss = fill(eps, length(dat[:betas]))
        scatter!(p, epss, dat[:betas]; 
            label = "", color = :black, m = 8, alpha = 0.7
        )
    end
    mysavefig(p, "beta_vs_eps"; MODsym) 
end
