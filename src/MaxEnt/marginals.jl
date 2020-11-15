function vatp_marginal_pdf(pol, beta; n = Int(1e5))
    
    vatp0, vatp1 = vatp_global_min(pol), vatp_global_max(pol)
    
    rvatp = range(vatp0, vatp1; length = n)
    Δvatp = step(rvatp)
    chuncks = zeros(n)
    
    f(vatp) = exp(-beta*((vatp1 - vatp)/pol.Nb))
    for (i, vatp) in enumerate(rvatp[1:end-1])
        chuncks[i] = f(vatp) * Δvg(vatp, pol)
    end
    Z = sum(chuncks)
    probs = chuncks ./ Z
    return rvatp, probs ./ Δvatp
end

function vg_marginal_pdf(pol, beta; n = Int(1e5), n2 = Int(1e2))
    
    vg0, vg1 = vg_global_min(pol), vg_global_max(pol)
    
    rvg = range(vg0, vg1; length = n)
    Δvg = step(rvg)
    chuncks = zeros(n)
    
    f(vatp) = exp(-beta*((vatp_global_max(pol) - vatp)/pol.Nb))
    @inbounds for (i, vg) in enumerate(rvg[1:end-1])
        vatp0, vatp1 = vatp_local_min(vg, pol), vatp_local_max(vg, pol)
        rvatp = range(vatp0, vatp1; length = n2)
        Δvatp = step(rvatp)
        for vatp in rvatp
            chuncks[i] += f(vatp)
        end
        chuncks[i] *= Δvatp
    end
    Z = sum(chuncks)
    probs = chuncks ./ Z
    return rvg, probs ./ Δvg
end

