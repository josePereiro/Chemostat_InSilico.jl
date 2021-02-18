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
dat_file(;sim_params...) = joinpath(
    InCh.DYN_DATA_DIR, 
    InLP.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...)
)
idxdat(dk, indexks...) = InLP.idxdat(DINDEX, dk, indexks...)


## ----------------------------------------------------------------------------
# Prepare network
const WLOCK = ReentrantLock()

const ME_HOMO = :ME_HOMO           # Do not use extra constraints
const ME_EXPECTED = :ME_EXPECTED   # Match ME and Dy biom average
const ME_CUTTED = :ME_CUTTED     # Match ME and Dy biom average and constraint av_ug
const ME_BOUNDED = :ME_BOUNDED     # Fix biom around observed

const FBA_OPEN = :FBA_OPEN
const FBA_BOUNDED = :FBA_BOUNDED

## ----------------------------------------------------------------------------
function run_ME!(M, MEmode; LP_cache, δ, δμ, DyBiom, verbose = true)
    
    # MaxEnt fun
    y = InLP.Y # atp/biomass yield
    maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)
    
    # M.net.ub[7] = 0.0 # ATPM PROBLEM Delete

    if MEmode == ME_CUTTED
        # Fix av_ug
        net = M.net
        net.ub[M.vg_idx] = min(M.Vg, M.cg * M.D/ M.X)
        net.ub[M.vl_idx] = min(M.Vl, M.cl * M.D/ M.X)
        L, U = InLP.fva(M.net)
        net.lb .= L; net.ub .= U
    end

    beta0 = 0.0
    if MEmode == ME_EXPECTED || MEmode == ME_CUTTED
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

            show_info = it == 1 || rem(it, 50) == 0 || 
                it == maxiters || abs(target - f)/target < th
            show_info && lock(WLOCK) do
                thid = threadid()
                @info("Grad Descent ", 
                    it, MEmode, 
                    target, f, 
                    (target - f), 
                    beta, thid
                ); println()
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

    # M.net.ub[7] = 0.0 # ATPM PROBLEM Delete

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
    min_vg = Inf
    # (max_vatp, min_vg) will maximize the yield
    for (vatpi, vatp) in vatp_range |> enumerate
        vg_range = vg_ranges[vatpi]
        isempty(vg_range) && continue
        if vatp > max_vatp
            max_vatp = vatp
            min_vg = minimum(vg_range)
        end
    end
    @assert !isinf(max_vatp)

    fbaf(vatp, vg) = (vatp == max_vatp && vg == min_vg) ? 1.0 : 0.0
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
# TODO: make script args
REDO_MAXENT = false
REDO_FBA = false

# ## ----------------------------------------------------------------------------
# # Compat
# let
#     δ = 0.08 # marginal discretization factor
#     δμ = 0.01 # ME_BOUNDED biomass variance
#     Vls, Ds, ϵs, τs = DINDEX[[:Vls, :Ds, :ϵs, :τs]]
#     params = Iterators.product(Vls, Ds, ϵs, τs)
#     for (Vl, D, ϵ, τ) in params
#         cfile = dat_file(;Vl, D, ϵ, τ, δ, δμ)
#         if isfile(cfile) # Check caches
#             @info("Doing", cfile)
#             MDAT = deserialize(cfile)
#             MDAT[:M0] = idxdat([:M], Vl, D, ϵ, τ)
#             serialize(cfile, MDAT)
#         end
#     end

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
    
    params = Iterators.product(Vls, Ds, ϵs, τs)
    Ch = Channel(1) do Ch_
        for (Vl, D, ϵ, τ) in params
            put!(Ch_, (Vl, D, ϵ, τ))
        end
    end
    N = length(params)
    
    @threads for thid in 1:nthreads()
        for (Vl, D, ϵ, τ) in Ch
            cfile = dat_file(;Vl, D, ϵ, τ, δ, δμ)
            
            ## ----------------------------------------------------------------------------
            MDAT = UJL.DictTree()
            M0 = MDAT[:M0] = idxdat([:M], Vl, D, ϵ, τ)
            status = idxdat([:status], Vl, D, ϵ, τ)
            LP_cache = nothing
            
            c = nothing
            lock(WLOCK) do
                gc += 1; c = gc
                LP_cache = InLP.vgvatp_cache(M0)
                push!(INDEX[:Vls], Vl); push!(INDEX[:Ds], D)
                push!(INDEX[:ϵs], ϵ); push!(INDEX[:τs], τ)
                INDEX[:DFILE, Vl, D, ϵ, τ] = relpath(cfile, InCh.PROJECT_DIR)
                INDEX[:STATUS, Vl, D, ϵ, τ] = status
                @info("Doing $c, prog: $gc/$N ... ", 
                    (Vl, D, ϵ, τ), 
                    M0.X, status, 
                    thid
                ); println()
            end
            
            ## ----------------------------------------------------------------------------
            if status != :stst # Only accept steady states
                lock(WLOCK) do
                    @info("Not a Stst (Skipping) $c, prog: $gc/$N ... ",
                        (Vl, D, ϵ, τ),
                        M0.X, status,
                        thid
                    ); println()
                end
                continue 
            end
            
            ## ----------------------------------------------------------------------------
            if isfile(cfile) # Check caches
                if REDO_MAXENT || REDO_FBA
                    MDAT = deserialize(cfile)
                    lock(WLOCK) do
                        @info("Cache found (Redoing) $c, prog: $gc/$N ... ", 
                            (Vl, D, ϵ, τ), 
                            status, basename(cfile),
                            REDO_MAXENT, REDO_FBA,
                            thid
                        ); println()
                    end
                else
                    lock(WLOCK) do
                        @info("Cache found (Skipping) $c, prog: $gc/$N ... ", 
                            (Vl, D, ϵ, τ), 
                            status, basename(cfile),
                            thid
                        ); println()
                    end
                    continue # skip
                end
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
            for MEmode in [ME_HOMO, ME_BOUNDED, ME_CUTTED, ME_EXPECTED]
                isfile(cfile) && !REDO_MAXENT && break

                    # Setup network
                    M = deepcopy(M0)
                    
                    MEMs, beta0 = run_ME!(M, MEmode; 
                        LP_cache, δ, δμ, DyBiom, verbose = false)
                    MEBiom = InLP.av(MEMs[InLP.BIOMASS_IDER]) # biomass dynamic mean

                    # Ranges
                    vatp_range, vg_ranges = InLP.vatpvg_ranges(M)

                    lock(WLOCK) do
                        MDAT[MEmode, :M] = M
                        MDAT[MEmode, :Ms] = MEMs
                        MDAT[MEmode, :beta0] = beta0
                        MDAT[MEmode, :POL] = (;vatp_range, vg_ranges)

                        @info("Done MaxEnt  $c, prog: $gc/$N ... ",
                            MEmode,  
                            (Vl, D, ϵ, τ),
                            M0.X, DyBiom, MEBiom,
                            thid
                        ); println()
                    end
            end # for MEmode

            ## ----------------------------------------------------------------------------
            # FBA
            for FBAmode in [FBA_BOUNDED, FBA_OPEN]
                isfile(cfile) && !REDO_FBA && break

                    # Setup network
                    M = deepcopy(M0)
                    
                    FBAMs = run_FBA!(M, FBAmode; 
                        LP_cache, δ, δμ, DyBiom, verbose = false)
                    FBABiom = InLP.av(FBAMs[InLP.BIOMASS_IDER]) # biomass dynamic mean

                    # Ranges
                    vatp_range, vg_ranges = InLP.vatpvg_ranges(M)

                    lock(WLOCK) do
                        MDAT[FBAmode, :M] = M
                        MDAT[FBAmode, :Ms] = FBAMs
                        MDAT[FBAmode, :POL] = (;vatp_range, vg_ranges)

                        @info("Done FBA  $c, prog: $gc/$N ... ",
                        FBAmode,  
                            (Vl, D, ϵ, τ),
                            M0.X, DyBiom, FBABiom,
                            thid
                        ); println()
                    end
            end

            ## ----------------------------------------------------------------------------
            lock(WLOCK) do
                @info("Finished  $c, prog: $gc/$N ... ",
                    (Vl, D, ϵ, τ),
                    M0.X, basename(cfile),
                    thid
                ); println()
                serialize(cfile, MDAT)
            end
            GC.gc()
        end # for (Vl, D, ϵ, τ) in Ch
    end #  @threads

    sort!(unique!(INDEX[:Vls])); sort!(unique!(INDEX[:Ds]))
    sort!(unique!(INDEX[:ϵs])); sort!(unique!(INDEX[:τs]))
end

## ----------------------------------------------------------------------------
# SAVING
UJL.save_data(InCh.MARGINALS_INDEX_FILE, INDEX)


