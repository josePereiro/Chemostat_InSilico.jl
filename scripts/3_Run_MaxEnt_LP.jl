import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Dates
    using Serialization
    using Random
end

## ----------------------------------------------------------------------------
# Load and clear DAT
# DINDEX [Vl, D, ϵ, τ] 
DINDEX = UJL.load_data(InCh.DYN_DATA_INDEX_FILE) # Dynamic index
DATA_FILE_PREFFIX = "marginal_dat"
dat_file(;sim_params...) = joinpath(InCh.DYN_DATA_DIR, 
    InLP.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...))
idxdat(dk, indexks...) = InLP.idxdat(DINDEX, dk, indexks...)

## ----------------------------------------------------------------------------
let
    DINDEX[:Ds]
end
## ----------------------------------------------------------------------------
# Prepare network
const WLOCK = ReentrantLock()
const STST_POL = :STST_POL   # Use polytope computed from chemostat stst assertion
const DYN_POL = :DYN_POL     # Use dynamic polytope

const ME_HOMO = :ME_HOMO           # Do not use extra constraints
const ME_EXPECTED = :ME_EXPECTED   # Match ME and Dy biom average
const ME_BOUNDED = :ME_BOUNDED     # Fix biom around observed
const FBA_OPEN = :FBA_OPEN
const FBA_BOUNDED = :FBA_BOUNDED

function prepare_pol!(M, POLTsym)
    net = M.net
    xi = M.X / M.D
    if POLTsym == DYN_POL
        # Use Michaelis-Menden
        net.ub[M.vg_idx] = max(net.lb[M.vg_idx], 
            min(M.cg / xi, (M.Vg * M.sg) / (M.Kg + M.sg)))
        net.ub[M.vl_idx] = max(net.lb[M.vl_idx], 
            min(M.cl / xi, (M.Vl * M.sl) / (M.Kl + M.sl)))
    elseif  POLTsym == STST_POL
        # Use Chemostat Steady-State assumption
        # ub < max(V, c/ xi) see cossio's paper
        net.ub[M.vg_idx] = max(net.lb[M.vg_idx], min(M.Vg , M.cg / xi))
        net.ub[M.vl_idx] = max(net.lb[M.vl_idx], min(M.Vl , M.cl / xi))
    else
        error("Unknown pol type")
    end
    return net
end

## ----------------------------------------------------------------------------
function run_ME!(M, MEmode; LP_cache, δ, δμ, DyBiom, verbose = true)
    
    # MaxEnt fun
    y = InLP.Y # atp/biomass yield
    maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)

    beta0 = 0.0
    if MEmode == ME_EXPECTED
        # Gradient descent
        target = DyBiom
        x0 = 1.5e3
        x1 = x0 * 0.9
        maxΔ = 100.0
        th = 1e-3
        maxiters = 500
        it = 1

        # Dynamic caching
        beta0 = UJL.grad_desc(;target, x0, x1, maxΔ, th, maxiters, verbose) do beta
            MEMs = InLP.get_marginals(maxentf(beta), M, [InLP.BIOMASS_IDER]; 
                δ, verbose = false, LP_cache)
            f = InLP.av(MEMs[InLP.BIOMASS_IDER])

            show_info = it == 1 || rem(it, 50) == 0 || it == maxiters
            if show_info
                lock(WLOCK) do
                    @info("Grad Descent ", 
                        it, target, 
                        f, (target - f), 
                        beta, threadid()
                    ); println()
                end
            end

            it += 1
            return f
        end
    end

    if MEmode == ME_BOUNDED
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = DyBiom * (1.0 + δμ)
        net.lb[M.obj_idx] = DyBiom * (1.0 - δμ)
        L, U = InLP.fva(M.net)
        net.lb .= L; net.ub .= U
    end

    MEMs = InLP.get_marginals(maxentf(beta0), M; δ, LP_cache, verbose)
    return MEMs, beta0
end

## ----------------------------------------------------------------------------
function run_FBA!(M, FBAmode; LP_cache, δ, δμ, DyBiom, verbose = true)

    if FBAmode == FBA_BOUNDED
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = DyBiom * (1.0 + δμ)
        net.lb[M.obj_idx] = DyBiom * (1.0 - δμ)
        L, U = InLP.fva(M.net)
        net.lb .= L; net.ub .= U
    end

    # Find maximum feasible vatp
    vatp_range, vg_ranges = InLP.vatpvg_ranges(M)
    max_vatp = -Inf
    for (vatpi, vatp) in vatp_range |> enumerate
        isempty(vg_ranges[vatpi]) && continue
        vatp > max_vatp && (max_vatp = vatp)
    end
    @assert !isinf(max_vatp)

    fbaf(vatp, vg) = vatp == max_vatp ? 1.0 : 0.0
    FBAMs = InLP.get_marginals(fbaf, M; δ, LP_cache, verbose)
    return FBAMs
end

# ## ----------------------------------------------------------------------------
# # Dev
# let
#     δ = 0.08
#     δμ = 0.01
#     Vl, D, ϵ, τ = 0.0, 0.036, 0.01, 0.0
#     FBAmode = :FBA_OPEN
#     POLTsym = :STST_POL
#     M = idxdat([:M], Vl, D, ϵ, τ)
#     prepare_pol!(M, POLTsym)

#     LP_cache = InLP.vgvatp_cache(M; marginf = 1.5)
#     # Dynamic marginal
#     f(vatp, vg) = M.Xb[vatp][vg] / M.X
#     DyMs = InLP.get_marginals(f, M; 
#         δ, LP_cache, verbose = false)
#     DyBiom = InLP.av(DyMs[InLP.BIOMASS_IDER]) # biomass dynamic mean
#     # @show DyBiom
    
#     # fba_sol = InLP.fba(M.net)
#     # vatp_fba = fba_sol[M.vatp_idx]
#     # vg_fba = fba_sol[M.vg_idx]
    
#     # # Discretization RoundDown ensure the values are inside the polytope
#     # vatp_fba = InLP.discretize(vatp_fba, M.δvatp; mode = RoundDown)
#     # vg_fba = InLP.discretize(vg_fba, M.δvg; mode = RoundDown)
#     # fbaf(vatp, vg) = (vatp == vatp_fba && vg == vg_fba) ? 1.0 : 0.0
#     # FBAMs = InLP.get_marginals(fbaf, M; δ, LP_cache, verbose = true)

#     run_FBA!(M, FBAmode; LP_cache, δ, δμ, DyBiom, verbose = false)

# end

## ----------------------------------------------------------------------------
# COMPUTE MARGINALS
INDEX = UJL.DictTree() # marginals dat
let
    δ = INDEX[:δ]   = 0.08 # marginal discretization factor
    δμ = INDEX[:δμ] = 0.01 # ME_BOUNDED biomass variance
    gc = 0

    INDEX[:Vls], INDEX[:Ds] = [], []
    INDEX[:ϵs], INDEX[:τs] = [], []
    
    Vls, Ds, ϵs, τs = DINDEX[[:Vls, :Ds, :ϵs, :τs]]
    # Ds = Ds[1:2:end] # Test
    params = Iterators.product(Vls, Ds, ϵs, τs) |> collect |> shuffle!
    N = length(params)
    @threads for (Vl, D, ϵ, τ) in params
        cfile = dat_file(;Vl, D, ϵ, τ, δ, δμ)
        
        ## ----------------------------------------------------------------------------
        MDAT = UJL.DictTree()
        M0 = idxdat([:M], Vl, D, ϵ, τ)
        status = idxdat([:status], Vl, D, ϵ, τ)
        LP_cache = nothing
        
        c = nothing
        lock(WLOCK) do
            gc += 1; c = gc
            LP_cache = InLP.vgvatp_cache(M0; marginf = 1.5)
            push!(INDEX[:Vls], Vl); push!(INDEX[:Ds], D)
            push!(INDEX[:ϵs], ϵ); push!(INDEX[:τs], τ)
            INDEX[:DFILE, Vl, D, ϵ, τ] = cfile
            INDEX[:STATUS, Vl, D, ϵ, τ] = status
            @info("Doing $c/$N ... ", 
                (Vl, D, ϵ, τ), 
                M0.X, status, 
                threadid()
            ); println()
        end
        
        ## ----------------------------------------------------------------------------
        if status != :stst # Only accept steady states
            lock(WLOCK) do
                @info("Not a Stst (Skipping) $c/$N ... ",
                    (Vl, D, ϵ, τ),
                    M0.X, status,
                    threadid()
                ); println()
            end
            continue 
        end
        
        ## ----------------------------------------------------------------------------
        if isfile(cfile) # Check caches
            lock(WLOCK) do
                @info("Cache found (Skipping) $c/$N ... ", 
                    (Vl, D, ϵ, τ), 
                    status, basename(cfile),
                    threadid()
                ); println()
            end
            continue
        end
        
        ## ----------------------------------------------------------------------------
        # Dynamic marginal
        f(vatp, vg) = M0.Xb[vatp][vg] / M0.X
        DyMs = InLP.get_marginals(f, M0; 
            δ, LP_cache, verbose = false)
        DyBiom = InLP.av(DyMs[InLP.BIOMASS_IDER]) # biomass dynamic mean
        lock(WLOCK) do
            MDAT[:DyMs] = DyMs
        end

        ## ----------------------------------------------------------------------------
        # MaxEnt marginals
        for MEmode in [ME_HOMO, ME_BOUNDED, ME_EXPECTED]
            for POLTsym in [STST_POL, DYN_POL]
                
                lock(WLOCK) do
                    @info("Doing  $c/$N ... ",
                        MEmode, POLTsym, 
                        (Vl, D, ϵ, τ),
                        M0.X, 
                        threadid()
                    ); println()
                end

                # Setup network
                M = deepcopy(M0)
                prepare_pol!(M, POLTsym)
                
                MEMs, beta0 = run_ME!(M, MEmode; 
                    LP_cache, δ, δμ, DyBiom, verbose = false)
                
                # Ranges
                vatp_range, vg_ranges = InLP.vatpvg_ranges(M)

                lock(WLOCK) do
                    MDAT[MEmode, POLTsym, :M] = M
                    MDAT[MEmode, POLTsym, :Ms] = MEMs
                    MDAT[MEmode, POLTsym, :beta0] = beta0
                    MDAT[MEmode, POLTsym, :POL] = (;vatp_range, vg_ranges)
                end
            end
        end

        ## ----------------------------------------------------------------------------
        # FBA
        for FBAmode in [FBA_BOUNDED, FBA_OPEN]
            for POLTsym in [STST_POL, DYN_POL]

                lock(WLOCK) do
                    @info("Doing  $c/$N ... ", 
                        FBAmode, POLTsym,
                        (Vl, D, ϵ, τ),
                        M0.X, 
                        threadid()
                    ); println()
                end

                # Setup network
                M = deepcopy(M0)
                prepare_pol!(M, POLTsym)
                
                FBAMs = run_FBA!(M, FBAmode; 
                    LP_cache, δ, δμ, DyBiom, verbose = false)
                
                # Ranges
                vatp_range, vg_ranges = InLP.vatpvg_ranges(M)

                lock(WLOCK) do
                    MDAT[FBAmode, POLTsym, :M] = M
                    MDAT[FBAmode, POLTsym, :Ms] = FBAMs
                    MDAT[FBAmode, POLTsym, :POL] = (;vatp_range, vg_ranges)
                end
            end
        end

        ## ----------------------------------------------------------------------------
        lock(WLOCK) do
            @info("Finished  $c/$N ... ",
                (Vl, D, ϵ, τ),
                M0.X, basename(cfile),
                threadid()
            ); println()
            serialize(cfile, MDAT)
        end
        GC.gc()
    end #  for (Vl, D, ϵ, τ)

    sort!(unique!(INDEX[:Vls])); sort!(unique!(INDEX[:Ds]))
    sort!(unique!(INDEX[:ϵs])); sort!(unique!(INDEX[:τs]))
end

## ----------------------------------------------------------------------------
# SAVING
UJL.save_data(InCh.MARGINALS_INDEX_FILE, INDEX)


