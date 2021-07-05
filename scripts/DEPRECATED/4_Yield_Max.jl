import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement
    const InU = InCh.Utilities

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Dates, Serialization, Random
    import JuMP, GLPK

    using Plots, GR
    GR.inline("png")
end

## ----------------------------------------------------------------------------
# Prepare network
const STST_POL = :STST_POL   # Use polytope computed from chemostat stst assertion
const DYN_POL = :DYN_POL     # Use dynamic polytope
const HOMO = :HOMO           # Do not use extra constraints
const EXPECTED = :EXPECTED         # Match ME and Dy biom average
const BOUNDED = :BOUNDED       # Fix biom around observed

MOD_COLORS = Dict(
    HOMO => :red, 
    BOUNDED => :orange,
    EXPECTED => :blue,
)

POL_STYLE = Dict(
    STST_POL => :dot,
    DYN_POL => :dash,
)

## ----------------------------------------------------------------------------
# Marginals
# MDAT[MODsym, POLTsym, :M, Vl, D, ϵ, τ]
# MDAT[MODsym, POLTsym, :MEMs, Vl, D, ϵ, τ]
# MDAT[MODsym, POLTsym, :beta0, Vl, D, ϵ, τ]
# MDAT[:STATUS, Vl, D, ϵ]
MINDEX = UJL.load_data(InCh.MARGINALS_INDEX_FILE)
POLTsym = STST_POL
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...);
idxdat(dk, indexks...; cache = true) = InLP.idxdat(MINDEX, dk, indexks...; cache)
WLOCK = ReentrantLock();
const fileid = "4"
mysavefig(p, pname; params...) = 
    InLP.mysavefig(p, pname, InLP.DYN_FIGURES_DIR, fileid; params...)

## ----------------------------------------------------------------------------
function avflx(M, idxs = eachindex(M.net.rxns); 
        LP_cache = nothing
    )

    if isnothing(LP_cache) 
        LP_cache = InLP.vgvatp_cache(M)
    end

    flxs = Dict(idx => 0.0 for idx in idxs)
    for (vatp, Xl) in M.Xb
        lcache = LP_cache[vatp]
        for (vg, X) in Xl
            p = X / M.X
            for idx in idxs
                flxs[idx] += lcache[vg][idx] * (X / M.X)
            end
        end
    end
    [flxs[idx] for idx in idxs]
end

## ----------------------------------------------------------------------------
ydat_file = joinpath(InLP.DYN_DATA_DIR, "ydat.jls")
YDAT = isfile(ydat_file) ? deserialize(ydat_file) : UJL.DictTree()
@time let
    # Vl = MINDEX[:Vls] |> first
    # τ = MINDEX[:τs] |> first
    # D = (MINDEX[:Ds] |> sort)[6]
    # ϵ = MINDEX[:ϵs][2]
    MODsym = HOMO
    POLTsym = STST_POL

    L = length(EXP_PARAMS)
    for (i, (Vl, D, ϵ, τ)) in EXP_PARAMS |> enumerate
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        haskey(YDAT, Vl, D, ϵ, τ) && continue # cached

        # prepare model
        M0 = idxdat([MODsym, POLTsym, :M], Vl, D, ϵ, τ; cache = false)
        net = M0.net |> deepcopy
        glcidx = InLP.rxnindex(net, "gt")
        objidx = InLP.rxnindex(net, "biom")
        M, N = size(net.S)

        # bound biomass
        dynflxs = avflx(M0)
        z = dynflxs[objidx]
        net.lb[objidx] = net.ub[objidx] = z

        # yield maximization
        yd = zeros(N); yd[glcidx] = 1.0
        net.c .= 0.0; net.c[objidx] = 1.0
        status, yflxs, yield = InLP.yLP(net, yd)

        @info("Yield max: ", (Vl, D, ϵ, τ), 
            z, status, yield, 
            dynflxs[objidx], yflxs[objidx], 
            threadid()
        ); println()
        YDAT[Vl, D, ϵ, τ] = (;net, status, dynflxs, yflxs, yield)

        up = i == L || rem(i, 5) == 0 
        up && serialize(ydat_file, YDAT) # caching
    end
end

## ----------------------------------------------------------------------------
# Plots
let
    expparams = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    ϵs = MINDEX[:ϵs]
    colors = Plots.distinguishable_colors(length(ϵs))
    rxns = ["gt", "vatp"]
    ps1 = Plots.Plot[]
    for (ϵi, ϵ) in ϵs |> enumerate
        ps2 = Plots.Plot[]
        for (rxni, rxn) in rxns |> enumerate
            m, M = Inf, -Inf
            p = Plots.plot(;title = string("rxn: ", rxn, " ϵ ", ϵ), xlabel = "dyn ave flxs", ylabel = "max yield flxs")
            for (i, (Vl, D, τ)) in expparams |> enumerate
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                net, status, dynflxs, yflxs, yield = YDAT[Vl, D, ϵ, τ]
                m = minimum([m, dynflxs[rxni], yflxs[rxni]])
                M = maximum([M, dynflxs[rxni], yflxs[rxni]])
                dm = abs(M - m) < 0.05 ? 0.05 : abs(M - m) * 0.1
                Plots.scatter!(p, [dynflxs[rxni]], [yflxs[rxni]]; 
                    label = "", color = colors[ϵi], ms = 7, 
                    alpha = 0.8,
                    xlim = [m - dm, M + dm],
                    ylim = [m - dm, M + dm],
                )
            end
            push!(ps2, p)
            push!(ps1, p)
        end
        mysavefig(ps2, "dyn_yield_tot_corrs"; ϵ)
    end
    layout = (4, length(rxns))
    mysavefig(ps1, "dyn_yield_tot_corrs"; layout)
end