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
                    ME_Z_OPEN_G_OPEN, 
                    # ME_Z_OPEN_G_BOUNDED, 
                    # ME_Z_EXPECTED_G_OPEN, 
                    ME_FULL_POLYTOPE,
                    ME_Z_EXPECTED_G_BOUNDED,
                    # ME_Z_EXPECTED_G_MOVING,
                    # ME_Z_FIXXED_G_OPEN, 
                    ME_Z_FIXXED_G_BOUNDED, 
                    ME_Z_EXPECTED_G_EXPECTED
                ]

                # check cache
                skip = REDO_MAXENT || !haskey(MDAT, MEmode)
                if skip; lock(WLOCK) do
                        @info("Cache found (Skipping) $c, prog: $gc/$N ... ", 
                            (Vl, D, ϵ, τ), MEmode, status, thid
                        ); println()
                    end; continue                
                end

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
                
                # check cache
                skip = REDO_FBA || !haskey(MDAT, FBAmode)
                if skip; lock(WLOCK) do
                        @info("Cache found (Skipping) $c, prog: $gc/$N ... ", 
                            (Vl, D, ϵ, τ), FBAmode, status, thid
                        ); println()
                    end; continue                
                end         

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
                DRY_RUN || serialize(cfile, MDAT) 
            end
            GC.gc()
        end # for (Vl, D, ϵ, τ) in Ch
    end #  @threads

    sort!(unique!(MINDEX[:Vls])); sort!(unique!(MINDEX[:Ds]))
    sort!(unique!(MINDEX[:ϵs])); sort!(unique!(MINDEX[:τs]))
end
