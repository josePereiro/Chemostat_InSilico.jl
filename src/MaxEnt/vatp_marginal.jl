# Marginals
function vatp_marginal(pol, beta; n = Int(1e5))
    
    vatp0, vatp1 = vatp_global_min(pol), vatp_global_max(pol)
    
    rvatp = range(vatp0, vatp1; length = n)
    Δvatp = step(rvatp)
    chuncks = zeros(n)
    
    f(vatp) = exp(-beta*((vatp1 - vatp)/pol.Nb))
    for (i, vatp) in enumerate(rvatp[1:end-1])
        chuncks[i] = f(vatp) * Δvg(vatp, pol) * Δvatp
    end
    Z = sum(chuncks)
    probs = chuncks ./ Z
    return rvatp, probs
end
