import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.MonteCarlo
using Plots
using Base.Threads
using ProgressMeter
using Test

## ------------------------------------------------------------------
# Marginals
function vatp_marginal_pdf_probs(pol, beta; n = Int(1e5))
    
    vatp0, vatp1 = vatpL(pol), vatpU(pol)
    
    rvatp = range(vatp0, vatp1; length = n)
    Δvatp = step(rvatp)
    chuncks = zeros(n)
    
    f(vatp) = exp(-beta*((vatp1 - vatp)/pol.Nb))
    for (i, vatp) in enumerate(rvatp[1:end-1])
        chuncks[i] = f(vatp) * Δvg(vatp, pol) * Δvatp
    end
    Z = sum(chuncks)
    probs = chuncks ./ Z
    @assert isapprox(sum(probs), 1.0; atol = 1e-7)
    return rvatp, probs
end


