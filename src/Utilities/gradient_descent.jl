function grad_desc(f;
    target::Vector,
    x0::Vector, x1::Vector,
    C::Vector = abs.(x0 - x1), th = 1e-5,
    maxiters = 1000, 
    toshow::Vector = [],
    Err = (fᵢ) -> abs.(target .- fᵢ) ./ ifelse.(iszero.(target), 1.0, abs.(target))
)

# initializing
xᵢ₋₁, xᵢ = x0, x1
fᵢ₋₁ = f(xᵢ₋₁)
ϵᵢ₋₁ = Err(fᵢ₋₁)
sense = ones(length(target))

prog = ProgressThresh(th, "Grad desc: ")
for it in 1:maxiters
    
    fᵢ = f(xᵢ)
    ϵᵢ = Err(fᵢ)
    sense .*= -sign.(ϵᵢ .- ϵᵢ₋₁)
    Δx = sense .* C .* ϵᵢ
    xᵢ += Δx

    maxϵᵢ = maximum(ϵᵢ)
    maxϵᵢ < th && break

    xᵢ₋₁ =  xᵢ
    ϵᵢ₋₁ = ϵᵢ

    update!(prog, maxϵᵢ; showvalues = vcat(
            [
                ("it", it),
                ("maxϵᵢ", maxϵᵢ),
                ("ϵᵢ", ϵᵢ),
                ("sense", sense),
                ("xᵢ", xᵢ),
                ("t", target),
                ("fᵢ", fᵢ),
            ], toshow
        )
    )
end
finish!(prog)

return xᵢ
end