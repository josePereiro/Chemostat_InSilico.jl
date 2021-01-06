discretize(n, Δ) = round(n / Δ) * Δ

function get_marginals(f::Function, M::SimModel, rxns = M.net.rxns; verbose = true, δ = 0.06, LP_cache = nothing)
    isempty(M.Xb) && error("Xb is empty!!!")
    
    if isnothing(LP_cache)
        LP_cache = vgvatp_cache(M)
    end
    
    # Collecting
    verbose && (prog = Progress(3 * length(rxns)))
    raw_marginals = Dict()
    for rxn in rxns
        rxni = rxnindex(M.net, rxn)
        raw_marginal = get!(raw_marginals, rxn, Dict{Float64, Float64}())
        for (vatp, lX) in M.Xb
            lcache = LP_cache[vatp]
            for (vg, X) in lX
                flx = lcache[vg][rxni]
                get!(raw_marginal, flx, 0.0)
                raw_marginal[flx] += f(vatp, vg)
            end
        end

        verbose && next!(prog)
    end

    # rounding
    marginals = Dict()
    for (rxn, raw_marginal) in raw_marginals
        marginal = get!(marginals, rxn, Dict{Float64, Float64}())
        
        flxs = keys(raw_marginal)
        flxL, flxU = minimum(flxs), maximum(flxs)
        Δ = abs(flxL - flxU) * δ

        for (raw_flx, p) in raw_marginal
            dflx = discretize(raw_flx, Δ)
            get!(marginal, dflx, 0.0)
            marginal[dflx] += p
        end

        verbose && next!(prog)
    end

    # normalizing
    for (rxn, marginal) in marginals
        Z = sum(values(marginal))
        for (flx, p) in marginal
            marginal[flx] = p/Z
        end
        @assert isapprox(sum(values(marginal)), 1.0; atol = 1e-8)
        verbose && next!(prog)
    end

    verbose && finish!(prog)
    return marginals
end

marginal_av(marginal) = sum(p * v for (v, p) in marginal)
marginal_va(marginal) = (μ = marginal_av(marginal); sum(p * ((v - μ)^2) for (v, p) in marginal))