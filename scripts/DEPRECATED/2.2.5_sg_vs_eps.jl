## ----------------------------------------------------------------------------
# sg vs eps
let
    Ds =  INDEX[:Ds]
    Vl = INDEX[:Vls] |> first
    τ =  INDEX[:τs] |> first
    ϵs = INDEX[:ϵs]

    p = plot(;title = "dynamic stst", xlabel = "ϵ", ylabel = "sg")
    ps = Plots.Plot[]
    for D in Ds
        xs, ys = [], []
        for ϵ in ϵs |> sort |> reverse
            M = idxdat([:M], Vl, D, ϵ, τ)
            status = idxdat([:status], Vl, D, ϵ, τ)

            @info("Doing", (Vl, D, ϵ, τ), M.X, status)
            status != :stst && continue

            sg = M.sg
            push!(xs, ϵ); push!(ys, sg)
        end
        length(ys) < length(ϵs) - 1 && continue
        plot!(p, xs, ys; label = "", alpha = 0.5, lw = 3)
    end
    mysavefig(p, "eps_vs_sg") 
end