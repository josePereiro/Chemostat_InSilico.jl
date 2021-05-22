# TODO: make this the default discretization form
_discretize(n, Δ) = iszero(Δ) ? n : round(n / Δ) * Δ

function get_marginals(f::Function, M::SimModel, rxns = M.net.rxns; 
        verbose = true, δ = 0.06, LP_cache = nothing)
    isempty(M.Xb) && error("Xb is empty!!!")
    
    isnothing(LP_cache) && (LP_cache = vgvatp_cache(M))
    
    vatp_range, vg_ranges = vatpvg_ranges(M)

    # COLLECTING
    verbose && (prog = Progress(3 * length(rxns)))
    raw_marginals = Dict{String, Dict{Float64, Float64}}()
    for rxn in rxns
        rxni = rxnindex(M.net, rxn)
        raw_marginal = get!(raw_marginals, rxn, Dict{Float64, Float64}())
        for (vatpi, vatp) in enumerate(vatp_range)
            lcache = LP_cache[vatp]
            for vg in vg_ranges[vatpi]
                flx = lcache[vg][rxni]
                get!(raw_marginal, flx, 0.0)
                # call function
                raw_marginal[flx] += f(vatp, vg)
            end
        end

        verbose && next!(prog)
    end

    # ROUNDING
    marginals = Dict{String, Dict{Float64, Float64}}()
    for (rxn, raw_marginal) in raw_marginals
        marginal = get!(marginals, rxn, Dict{Float64, Float64}())
        
        flxs = keys(raw_marginal)
        flxL, flxU = minimum(flxs), maximum(flxs)
        Δ = abs(flxL - flxU) * δ
        for (raw_flx, p) in raw_marginal
            dflx = _discretize(raw_flx, Δ)
            get!(marginal, dflx, 0.0)
            marginal[dflx] += p
        end

        verbose && next!(prog)
    end

    # NORMALIZING
    for (_, marginal) in marginals
        Z = sum(values(marginal))
        @assert Z > 0
        for (flx, p) in marginal
            marginal[flx] = p/Z
        end

        # sum_ = sum(values(marginal))
        # @assert(isapprox(sum_, 1.0; atol = 1e-8), sum_)
        verbose && next!(prog)
    end

    verbose && finish!(prog)
    return marginals
end

function Δv(marginal)
    sum = 0
    for v in keys(marginal)
        sum += v
    end
    sum/length(marginal)
end

av(marginal) = sum(p * v for (v, p) in marginal) / Δv(marginal)
va(marginal) = (μ = av(marginal); sum(p * ((v - μ)^2) for (v, p) in marginal) / Δv(marginal))