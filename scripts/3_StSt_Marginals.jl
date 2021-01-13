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
    using Dates
    using Serialization
    using Random
end

## ----------------------------------------------------------------------------
let
    M0 = InLP.SimModel()
    box = InLP.pol_box(M0)
    @show box
    LP_cache = InLP.vgvatp_cache(M0; marginf = 1.5)
    vatps = keys(LP_cache)
    vatpbs = (minimum(vatps), maximum(vatps))
    @show vatpbs
end

## ----------------------------------------------------------------------------
# Load and clear DAT
# DINDEX [Vl, D, ϵ, τ] 
DINDEX = UJL.load_data(InCh.DYN_DATA_INDEX_FILE) # Dynamic index
DATA_FILE_PREFFIX = "marginal_dat"
dat_file(;sim_params...) = joinpath(InCh.DYN_DATA_DIR, 
    InLP.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...))

function getdat(dk, indexks...)
    FILE = DINDEX[:DFILE, indexks...]
    if FILE isa UJL.ITERABLE
        dat = []
        for F in FILE
            datum = deserialize(F)[dk]
            push!(dat, datum)
        end
        return dat
    else
        dat = deserialize(FILE)
        return dat[dk]
    end
end

## ----------------------------------------------------------------------------
# Prepare network
const WLOCK = ReentrantLock()
const STST_POL = :STST_POL   # Use polytope computed from chemostat stst assertion
const DYN_POL = :DYN_POL     # Use dynamic polytope
const HOMO = :HOMO           # Do not use extra constraints
const HETER = :HETER         # Match ME and Dy biom average
const FIXXED = :FIXXED       # Fix biom around observed

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

function run_model!(M, MODsym; LP_cache, δ, δμ, DyBiom, verbose = true)
    
    # MaxEnt fun
    y = InLP.Y # atp/biomass yield
    maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)

    if MODsym == HETER
        # Gradient descent
        target = DyBiom
        x0 = 1.5e3
        x1 = x0 * 0.9
        maxΔ = 100.0
        th = 1e-3
        maxiters = 500
        it = 1

        # Dynamic caching
        beta0 = InU.grad_desc(;target, x0, x1, maxΔ, th, maxiters, verbose) do beta
            MEMs = InLP.get_marginals(maxentf(beta), M, [InLP.BIOMASS_IDER]; 
                δ, verbose = false, LP_cache)
            f = InLP.av(MEMs[InLP.BIOMASS_IDER])
            if it == 1 || rem(it, 50) == 0 || it == maxiters
                lock(WLOCK) do
                    @info "Grad Descent " it target f (target - f) beta threadid()
                    println()
                end
            end
            it += 1
            return f
        end
    else
        beta0 = 0.0
    end

    if MODsym == FIXXED
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = DyBiom * (1.0 + δμ)
        net.lb[M.obj_idx] = DyBiom * (1.0 - δμ)
        L, U = InLP.fva(M.net)
        net.lb .= L
        net.ub .= U
    end

    MEMs = InLP.get_marginals(maxentf(beta0), M; δ, LP_cache, verbose)
    return MEMs, beta0
end

## ----------------------------------------------------------------------------
# COMPUTE MARGINALS
INDEX = UJL.DictTree() # marginals dat
let

    δ = INDEX[:δ]   = 0.08 # marginal discretization factor
    δμ = INDEX[:δμ] = 0.01 # FIXXED biomass variance

    N = prod(length.(DINDEX[[:Vls, :Ds, :ϵs, :τs]]))
    c = 0

    INDEX[:Vls] = []
    INDEX[:Ds] = []
    INDEX[:ϵs] = []
    INDEX[:τs] = []
    
    params = Iterators.product(DINDEX[[:Vls, :Ds, :ϵs, :τs]]...) |> collect |> shuffle!
    @threads for (Vl, D, ϵ, τ) in params
        cfile = dat_file(;Vl, D, ϵ, τ, δ, δμ)
        
        ## ----------------------------------------------------------------------------
        MDAT = UJL.DictTree()
        M0 = getdat(:M, Vl, D, ϵ, τ)
        LP_cache = nothing
        
        lock(WLOCK) do
            c += 1
            LP_cache = InLP.vgvatp_cache(M0; marginf = 1.5)
            push!(INDEX[:Vls], Vl)
            push!(INDEX[:Ds], D)
            push!(INDEX[:ϵs], ϵ)
            push!(INDEX[:τs], τ)
            INDEX[:DFILE, Vl, D, ϵ, τ] = cfile
            INDEX[:STATUS, Vl, D, ϵ, τ] = M0.X < DINDEX[:death_th] ? :death : :stst
            @info "Doing $c/$N ... " Vl, D, ϵ, τ M0.X threadid()
            println()
        end
        INDEX[:STATUS, Vl, D, ϵ, τ] == :death && continue  # Exclude deaths
        
        ## ----------------------------------------------------------------------------
        if isfile(cfile)
            lock(WLOCK) do
                @info "Cache found (Skipping)" Vl, D, ϵ, τ basename(cfile) threadid()
                println()
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
        for MODsym in [HOMO, FIXXED, HETER]
            for POLTsym in [STST_POL, DYN_POL]
                
                lock(WLOCK) do
                    @info "Doing  $c/$N ... " MODsym POLTsym Vl, D, ϵ, τ M0.X threadid()
                    println()
                end

                # Setup network
                M = deepcopy(M0)
                prepare_pol!(M, POLTsym)
                
                MEMs, beta0 = run_model!(M, MODsym; 
                    LP_cache, δ, δμ, DyBiom, verbose = false)
                
                # Ranges
                vatp_range, vg_ranges = InLP.vatpvg_ranges(M)

                lock(WLOCK) do
                    MDAT[MODsym, POLTsym, :M] = M
                    MDAT[MODsym, POLTsym, :MEMs] = MEMs
                    MDAT[MODsym, POLTsym, :beta0] = beta0
                    MDAT[MODsym, POLTsym, :POL] = (;vatp_range, vg_ranges)
                end
            end
        end

        lock(WLOCK) do
            @info "Finished  $c/$N ... " Vl, D, ϵ, τ M0.X basename(cfile) threadid()
            serialize(cfile, MDAT)
            println()
        end
        GC.gc()

    end #  for (Vl, D, ϵ, τ)

    unique!(INDEX[:Vls])
    unique!(INDEX[:Ds])
    unique!(INDEX[:ϵs])
    unique!(INDEX[:τs])
end

## ----------------------------------------------------------------------------
# SAVING
UJL.save_data(InCh.MARGINALS_INDEX_FILE, INDEX)


