import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Dates
    using Serialization
    using Random
end

# Compat
# @eval InCh const LP_Implement = Dynamic

## ----------------------------------------------------------------------------
# Load and clear DAT
# DINDEX [Vl, D, ϵ, τ] 
DINDEX_FILE = Dyn.procdir("dym_dat_index.bson")
DINDEX = UJL.load_data(DINDEX_FILE) # Dynamic index
DATA_FILE_PREFFIX = "marginal_dat"
dat_file(;sim_params...) = 
    Dyn.procdir(UJL.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...)) 
idxdat(dk, indexks...; cache = false, emptycache = true) = 
    Dyn.idxdat(DINDEX, dk, indexks...; cache, emptycache)

# ----------------------------------------------------------------------------
const WLOCK = ReentrantLock()

# ----------------------------------------------------------------------------
# METHOD VARIANTS
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_OPEN_G_BOUNDED       = :ME_Z_OPEN_G_BOUNDED        # 

const ME_Z_EXPECTED_G_OPEN      = :ME_Z_EXPECTED_G_OPEN       # Match ME and Dy biom average
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED   # 
const ME_FULL_POLYTOPE          = :ME_FULL_POLYTOPE           # 
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug

const ME_Z_FIXXED_G_OPEN        = :ME_Z_FIXXED_G_OPEN         # Fix biom around observed
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

const FBA_Z_OPEN_G_OPEN       = :FBA_Z_OPEN_G_OPEN
const FBA_Z_OPEN_G_BOUNDED    = :FBA_Z_OPEN_G_BOUNDED
const FBA_Z_FIXXED_G_OPEN     = :FBA_Z_FIXXED_G_OPEN 
const FBA_Z_FIXXED_G_BOUNDED  = :FBA_Z_FIXXED_G_BOUNDED

## ----------------------------------------------------------------------------
function run_ME!(M, MEmode; LP_cache, δ, δμ, biom_avPX, vg_avPX)
    
    thid = threadid()

    ## -----------------------------------------------------------
    # Globals
    cgD_X = M.cg * M.D/ M.X
    biom_idx = M.obj_idx
    z(vatp, vg) = LP_cache[vatp][vg][biom_idx]
    vg(vatp, vg) = vg
    
    ## -----------------------------------------------------------
    # G BOUNDING

    ## -----------------------------------------------------------
    ismode = MEmode in [ME_Z_OPEN_G_BOUNDED, ME_Z_EXPECTED_G_BOUNDED, ME_Z_FIXXED_G_BOUNDED]
    ismode && let
        # Fix av_ug
        net = M.net
        net.ub[M.vg_idx] = min(M.Vg, cgD_X)
        net.ub[M.vl_idx] = min(M.Vl, cgD_X)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    ## -----------------------------------------------------------
    # Z FIXXED

    ## -----------------------------------------------------------
    ismode = MEmode in [ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED]
    ismode && let
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = biom_avPX * (1.0 + δμ)
        net.lb[M.obj_idx] = biom_avPX * (1.0 - δμ)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    ## -----------------------------------------------------------
    # EXPECTED
    
    ## -----------------------------------------------------------
    # Globals
    beta_biom = 0.0
    beta_vg = 0.0
    maxiter = 800
    verb_frec = 50.0
    gdth = 1e-2

    ## -----------------------------------------------------------
    ismode = MEmode in [ME_Z_EXPECTED_G_OPEN, ME_Z_EXPECTED_G_BOUNDED]
    ismode && let
        # Gradient descent
        target = biom_avPX
        x0 = 1.5e2
        x1 = x0 * 0.9
        maxΔx = 100.0
        gdit = 1

        # grad desc
        join = Dyn.get_join(M)
        function f1(gdmodel)
            
            beta_biom = UJL.gd_value(gdmodel)

            PME = Dyn.get_join!(M, join) do vatp, vg
                exp(beta_biom * z(vatp, vg))
            end
            biom_avPME = Dyn.ave_over(z, PME)
            
            biom_err = gdmodel.ϵi
            show_info = gdit == 1 || rem(gdit, verb_frec) == 0 || 
                gdit == maxiter || biom_err < gdth
            show_info && lock(WLOCK) do
                @info("Grad Descent ", 
                    gdit, MEmode, 
                    (biom_avPX, biom_avPME), 
                    biom_err, beta_biom, thid
                ); println()
            end

            gdit += 1
            return biom_avPME
        end

        gdmodel = UJL.grad_desc(f1; gdth, 
            target, x0, x1, maxΔx, maxiter, 
            verbose = false
        )
        beta_biom = UJL.gd_value(gdmodel)
    end

    ## -----------------------------------------------------------
    ismode = MEmode == ME_Z_EXPECTED_G_EXPECTED
    ismode && let
        target = [biom_avPX, vg_avPX]
        x0 = [100.0, -10.0]
        x1 = [101.0, -11.0]
        maxΔx = [80.0, 30.0]
        gdit = 1

        join = Dyn.get_join(M)
        function up_fun!(gdmodel)

            beta_biom, beta_vg = UJL.gd_value(gdmodel)
            PME = Dyn.get_join!(M, join) do vatp_, vg_
                exp(beta_biom * z(vatp_, vg_) + beta_vg * vg(vatp_, vg_))
            end
            biom_avPME = Dyn.ave_over(z, PME)
            vg_avPME = Dyn.ave_over(vg, PME)
            
            biom_err = abs(biom_avPX - biom_avPME)/biom_avPX
            vg_err = abs(vg_avPX - vg_avPME)/vg_avPX
            err = max(biom_err, vg_err)

            show_info = gdit == 1 || rem(gdit, verb_frec) == 0 || 
                gdit == maxiter || err < gdth
            show_info && lock(WLOCK) do
                @info("Grad Descent ", 
                    gdit, MEmode, 
                    (biom_avPX, biom_avPME),
                    (vg_avPX, vg_avPME),
                    err, beta_biom, beta_vg, 
                    thid
                ); println()
            end

            gdit += 1
            return [biom_avPME, vg_avPME]
        end

        gdmodel = UJL.grad_desc_vec(up_fun!; 
            target, x0, x1, maxΔx, gdth, maxiter, 
            verbose = false
        )
        beta_biom, beta_vg = UJL.gd_value(gdmodel)
    end

    ## -----------------------------------------------------------
    ismode = MEmode == ME_FULL_POLYTOPE
    ismode && let
        topiter = 1
        join = Dyn.get_join(M)
        stth, stw = 0.1, 8
        betas_biom, betas_vg = [beta_biom], [beta_vg]
        biom_avPME, vg_avPME = 0.0, 0.0

        while true

            ## -----------------------------------------------------------------------------------------------
            # find z_beta⁰
            # Gradient descent
            target = biom_avPX
            x0 = beta_biom
            maxΔx = max(100.0, abs(beta_biom) * 0.1)
            x1 = x0 + maxΔx * 0.1

            function z_fun(gdmodel)

                gditer = gdmodel.iter
                beta_biom = UJL.gd_value(gdmodel)

                PME = Dyn.get_join!(M, join) do vatp_, vg_
                    exp(beta_biom * z(vatp_, vg_) + beta_vg * vg(vatp_, vg_))
                end
                biom_avPME = Dyn.ave_over(z, PME)
                
                err = gdmodel.ϵi
                show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
                    gditer == maxiter || err < gdth
                show_info && begin
                    @info("z grad Descent ", 
                        topiter, gditer, 
                        MEmode, 
                        (biom_avPX, biom_avPME), 
                        err, beta_biom, thid
                    ); println()
                end

                return biom_avPME 
            end

            gdmodel = UJL.grad_desc(z_fun; gdth, 
                target, x0, x1, maxΔx, 
                maxiter, verbose = false
            )
            beta_biom = UJL.gd_value(gdmodel)
            
            ## -----------------------------------------------------------------------------------------------
            # Check balance at beta_vg = 0.0
            PME = Dyn.get_join!(M, join) do vatp_, vg_
                beta_vg_ = 0.0
                exp(beta_biom * z(vatp_, vg_) + beta_vg_ * vg(vatp_, vg_))
            end
            vg_avPME_at_b0 = Dyn.ave_over(PME) do vatp_, vg_
                vg(vatp_, vg_)
            end
            vg_valid_at_b0 = vg_avPME_at_b0 <= cgD_X

            ## -----------------------------------------------------------------------------------------------
            # if not valid, move beta_vg
            if vg_valid_at_b0
                beta_vg = 0.0
            else
                # Gradient descent
                target = cgD_X * (1.0 - gdth)
                x0 = beta_vg
                maxΔx = max(10.0, abs(beta_vg) * 0.1)
                x1 = x0 + maxΔx * 0.1
                
                function vg_fun(gdmodel)

                    gditer = gdmodel.iter
                    beta_vg = UJL.gd_value(gdmodel)

                    PME = Dyn.get_join!(M, PME) do vatp_, vg_
                        exp(beta_biom * z(vatp_, vg_) + beta_vg * vg(vatp_, vg_))
                    end
                    vg_avPME = Dyn.ave_over(vg, PME)
                    
                    err = gdmodel.ϵi
                    show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
                        gditer == maxiter || err < gdth
                    show_info && begin
                        @info("vg grad descent ", 
                            topiter, gditer,  
                            MEmode, 
                            (cgD_X, vg_avPME), 
                            err, beta_biom, thid
                        ); println()
                    end

                    return vg_avPME 
                end

                gdmodel = UJL.grad_desc(vg_fun; gdth, 
                    target, x0, x1, maxΔx, 
                    maxiter, verbose = false
                )
                beta_vg = UJL.gd_value(gdmodel)
            end

            push!(betas_biom, beta_biom); push!(betas_vg, beta_vg)

            conv = (vg_valid_at_b0 && gdmodel.ϵi < gdth) || 
                (UJL.is_stationary(betas_biom, stth, stw) && UJL.is_stationary(betas_vg, stth, stw))
            
            @info("Finished Round", 
                topiter, conv,
                MEmode,
                (beta_biom, beta_vg),
                (vg_avPME, cgD_X),
                (vg_avPME_at_b0, cgD_X),
                (biom_avPME, biom_avPX),
                thid
            ); println()

            conv && break
            
            topiter += 1
            topiter > maxiter && break
        end
    end

    ## -----------------------------------------------------------
    ismode = MEmode == ME_Z_EXPECTED_G_MOVING
    ismode && let

        # init globals
        Δstep = 0.5
        maxiter = 500
        PME = nothing
        biom_avPME = 0.0
        biom_err = Inf
        
        ## -----------------------------------------------------------
        last_beta = 50.0
        for rit in 1:maxiter

            ## -----------------------------------------------------------
            # Find biom beta
            target = biom_avPX
            x0 = beta_biom
            x1 = x0 * 0.9
            maxΔx = 100.0
        
            # grad desc
            it = 1
            join = Dyn.get_join(M)
            function up_fun!(gdmodel)
                beta = UJL.gd_value(gdmodel)

                PME = Dyn.get_join!(M, join) do vatp, vg
                    exp(beta * z(vatp, vg))
                end
                biom_avPME = Dyn.ave_over(z, PME)
                
                biom_err = gdmodel.ϵi
                show_info = it == 1 || rem(it, verb_frec) == 0 || 
                    it == maxiter || biom_err < gdth
                show_info && lock(WLOCK) do
                    @info("Grad Descent ", 
                        it, MEmode, 
                        (biom_avPX, biom_avPME), 
                        biom_err, beta, thid
                    ); println()
                end

                it += 1
                return biom_avPME
            end

            gdmodel = UJL.grad_desc(up_fun!; gdth, 
                target, x0, x1, maxΔx, maxiter, 
                verbose = false
            )
            beta_biom = UJL.gd_value(gdmodel)

            ## -----------------------------------------------------------
            # move vg
            vg_avPME = Dyn.ave_over(vg, PME)
        
            Δ = cgD_X - vg_avPME
            net = M.net
            net.ub[M.vg_idx] = max(net.lb[M.vg_idx], 
                min(M.Vg, net.ub[M.vg_idx] + Δstep * Δ)
            )
            L, U = Dyn.fva(M.net)
            net.lb .= L; net.ub .= U
            vg_ub = net.ub[M.vg_idx]
        
            valid_vg_avPME = Δ >= 0

            ## -----------------------------------------------------------
            lock(WLOCK) do 
                @info("End round", 
                    rit, beta_biom, 
                    biom_err, valid_vg_avPME, 
                    (M.Vg, vg_ub),
                    (cgD_X, vg_avPME),
                    thid
                ); println()
            end

            conv = valid_vg_avPME && biom_err < gdth
            conv &&break

        end # for rit in 1:maxiter
    end

    ## -----------------------------------------------------------
    # MARGINALS
    MEMs = Dyn.get_marginals(M; δ, LP_cache, verbose = false) do vatp_, vg_
        exp((beta_biom * z(vatp_, vg_)) + (beta_vg * vg(vatp_, vg_)))
    end
    return MEMs, beta_biom, beta_vg
end

## ----------------------------------------------------------------------------
function run_FBA!(M, FBAmode; LP_cache, δ, δμ, biom_avPX, verbose = true)

    # Z FIXXED
    ismode = FBAmode in [FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED]
    ismode && let
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = biom_avPX * (1.0 + δμ)
        net.lb[M.obj_idx] = biom_avPX * (1.0 - δμ)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    # G BOUNDING
    ismode = FBAmode in [FBA_Z_OPEN_G_BOUNDED, FBA_Z_FIXXED_G_BOUNDED]
    ismode && let
        # Fix av_ug
        net = M.net
        net.ub[M.vg_idx] = min(M.Vg, M.cg * M.D/ M.X)
        net.ub[M.vl_idx] = min(M.Vl, M.cl * M.D/ M.X)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    # FBA
    # Find maximum feasible vatp
    vatp_range, vg_ranges = Dyn.vatpvg_ranges(M)
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
    FBAMs = Dyn.get_marginals(fbaf, M; δ, LP_cache, verbose)
    return FBAMs
end

## ----------------------------------------------------------------------------
# TODO: make script args
REDO_MAXENT = false
REDO_FBA = false

## ----------------------------------------------------------------------------
# COMPUTE MARGINALS
MINDEX = UJL.DictTree() # marginals dat
let
    δ = MINDEX[:δ]   = 0.08 # marginal discretization factor
    δμ = MINDEX[:δμ] = 0.01 # ME_Z_FIXXED_G_OPEN biomass variance
    gc = 0

    MINDEX[:Vls], MINDEX[:Ds] = [], []
    MINDEX[:ϵs], MINDEX[:τs] = [], []
    
    Vls, Ds, ϵs, τs = DINDEX[[:Vls, :Ds, :ϵs, :τs]]
    
    params = Iterators.product(Vls, Ds, ϵs, τs)
    Ch = Channel(1) do Ch_
        for (Vl, D, ϵ, τ) in params
            put!(Ch_, (Vl, D, ϵ, τ))
        end
    end
    N = length(params)
    
    LP_cache = nothing
    @threads for thid in 1:nthreads()
        for (Vl, D, ϵ, τ) in Ch
            cfile = dat_file(;Vl, D, ϵ, τ, δ, δμ)
            
            ## ----------------------------------------------------------------------------
            MDAT = UJL.DictTree()
            M0 = MDAT[:M0] = idxdat([:M], Vl, D, ϵ, τ; cache = true)
            status = idxdat([:status], Vl, D, ϵ, τ; cache = false, emptycache = true)
            
            # LP_cache = nothing
            c = nothing
            lock(WLOCK) do
                gc += 1; c = gc
                
                isnothing(LP_cache) && (LP_cache = Dyn.vgvatp_cache(M0))

                push!(MINDEX[:Vls], Vl); push!(MINDEX[:Ds], D)
                push!(MINDEX[:ϵs], ϵ); push!(MINDEX[:τs], τ)
                MINDEX[:DFILE, Vl, D, ϵ, τ] = relpath(cfile, InCh.projectdir())
                MINDEX[:STATUS, Vl, D, ϵ, τ] = status
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
            fX(vatp, vg) = M0.Xb[vatp][vg] / M0.X
            DyMs = Dyn.get_marginals(fX, M0; δ, LP_cache, verbose = false)
            biom_avPX = Dyn.av(DyMs[Dyn.BIOMASS_IDER]) # biomass dynamic mean
            vg_avPX = Dyn.av(DyMs["gt"]) # biomass dynamic mean
            lock(WLOCK) do
                MDAT[:DyMs] = DyMs
            end
            DyMs = nothing

            ## ----------------------------------------------------------------------------  
            # MaxEnt marginals
            for MEmode in [ 
                    # ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
                    # ME_Z_EXPECTED_G_OPEN, 
                    ME_FULL_POLYTOPE,
                    # ME_Z_EXPECTED_G_BOUNDED,
                    # ME_Z_EXPECTED_G_MOVING,
                    # ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
                    # ME_Z_EXPECTED_G_EXPECTED
                ]

                isfile(cfile) && !REDO_MAXENT && break

                # Setup network
                M = deepcopy(M0)
                
                MEMs, beta_biom, beta_vg = run_ME!(M, MEmode; 
                    LP_cache, δ, δμ, 
                    biom_avPX, vg_avPX
                )
                biom_avPME = Dyn.av(MEMs[Dyn.BIOMASS_IDER]) # biomass dynamic mean
                vg_avPME = Dyn.av(MEMs["gt"]) # biomass dynamic mean

                # Ranges
                vatp_range, vg_ranges = Dyn.vatpvg_ranges(M)

                lock(WLOCK) do
                    MDAT[MEmode, :M] = M
                    MDAT[MEmode, :Ms] = MEMs
                    MDAT[MEmode, :beta_biom] = beta_biom
                    MDAT[MEmode, :beta_vg] = beta_vg
                    MDAT[MEmode, :POL] = (;vatp_range, vg_ranges)

                    @info("Done MaxEnt  $c, prog: $gc/$N ... ",
                        MEmode,  
                        (Vl, D, ϵ, τ),
                        M0.X, 
                        (biom_avPX, biom_avPME),
                        (vg_avPX, vg_avPME),
                        thid
                    ); println()
                end
            end # for MEmode

            ## ----------------------------------------------------------------------------
            # FBA
            for FBAmode in [
                    FBA_Z_OPEN_G_OPEN, FBA_Z_OPEN_G_BOUNDED,
                    FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED
                ]
                
                isfile(cfile) && !REDO_FBA && break

                # Setup network
                M = deepcopy(M0)
                
                FBAMs = run_FBA!(M, FBAmode; 
                    LP_cache, δ, δμ, biom_avPX, verbose = false)
                biom_avFBA = Dyn.av(FBAMs[Dyn.BIOMASS_IDER]) # biomass dynamic mean

                # Ranges
                vatp_range, vg_ranges = Dyn.vatpvg_ranges(M)

                lock(WLOCK) do
                    MDAT[FBAmode, :M] = M
                    MDAT[FBAmode, :Ms] = FBAMs
                    MDAT[FBAmode, :POL] = (;vatp_range, vg_ranges)

                    @info("Done FBA  $c, prog: $gc/$N ... ",
                    FBAmode,  
                        (Vl, D, ϵ, τ),
                        M0.X, 
                        (biom_avPX, biom_avFBA),
                        thid
                    ); println()
                end
            end

            ## ----------------------------------------------------------------------------
            # SAVING
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

    sort!(unique!(MINDEX[:Vls])); sort!(unique!(MINDEX[:Ds]))
    sort!(unique!(MINDEX[:ϵs])); sort!(unique!(MINDEX[:τs]))
end

## ----------------------------------------------------------------------------
# SAVING
MINDEX_FILE = Dyn.procdir("marg_dat_index.bson")
UJL.save_data(MINDEX_FILE, MINDEX)


