# ----------------------------------------------------------------------------
# beta vs eps
let
    method = ME_FULL_POLYTOPE
    ϵs = MINDEX[:ϵs] |> sort
    colors = Plots.distinguishable_colors(length(MINDEX[:Ds]))
    colors = Dict(D => c for (D, c) in zip(MINDEX[:Ds], colors))
    p = plot(;xlabel = "beta", ylabel = "ϵ")
    sparams = (;alpha = 0.5, ms = 6)
    exp_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    for (Vl, D, τ) in exp_params
        ϵ_ser = []
        beta_ser = []
        for ϵ in ϵs
            MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
            beta = idxdat([method, :beta_biom], Vl, D, ϵ, τ)
            push!(ϵ_ser, ϵ)
            push!(beta_ser, beta)
        end
        scatter!(p, beta_ser, ϵ_ser; label = "", color = colors[D], sparams...)
    end
    mysavefig(p, "biom_beta_vs_eps_D_colored"; method)
end
