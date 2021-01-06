function grad_desc(f;
        target::Float64,
        x0::Float64, x1::Float64,
        maxΔ::Float64 = abs(x0 - x1), th = 1e-5,
        maxiters = Int(1e4), 
        verbose = true,
        toshow::Vector = [],
        Err = (fᵢ) -> abs(target - fᵢ) / ifelse(iszero(target), 1.0, abs(target))
    )

    # initializing
    xᵢ₋₁, xᵢ = x0, x1
    fᵢ₋₁ = f(xᵢ₋₁)
    ϵᵢ₋₁ = Err(fᵢ₋₁)
    sense = one(target)

    verbose && (prog = ProgressThresh(th, "Grad desc: "))
    for it in 1:maxiters
        
        fᵢ = f(xᵢ)
        ϵᵢ = Err(fᵢ)
        sense *= ϵᵢ > ϵᵢ₋₁ ? -1.0 : 1.0
        Δx = sense * maxΔ * ϵᵢ
        xᵢ += Δx

        
        maxϵᵢ = maximum(ϵᵢ)
        maxϵᵢ < th && break

        (iszero(Δx) || isnan(Δx) || isinf(xᵢ)) && break

        xᵢ₋₁ =  xᵢ
        ϵᵢ₋₁ = ϵᵢ

        verbose && update!(prog, maxϵᵢ; showvalues = vcat(
                [
                    ("it", it),
                    ("maxϵᵢ", maxϵᵢ),
                    ("ϵᵢ", ϵᵢ),
                    ("sense", sense),
                    ("xᵢ", xᵢ),
                    ("Δx", Δx),
                    ("t", target),
                    ("fᵢ", fᵢ),
                ], toshow
            )
        )
    end
    verbose && finish!(prog)

    return xᵢ
end