let
    
    ME_MODELS = [
        ME_Z_OPEN_G_OPEN, 
        ME_FULL_POLYTOPE,
        ME_Z_EXPECTED_G_EXPECTED,
        ME_Z_EXPECTED_G_BOUNDED,
        ME_Z_FIXXED_G_BOUNDED, 
    ]
    
    # ---------------------------------------------------------------
    # COLLECT DATA
    # store collected data from join distributions
    nths = 3 # working threads
    JDAT_CID =("JOIN DATA CACHE")
    # UJL.delete_cache(JDAT_CID) # Reset cache
    JDAT = UJL.load_cache(JDAT_CID) do
        
        JDAT_ = Dict()
        LP_cache = nothing
        WLOCK = ReentrantLock()

        Ch = Channel(nthreads()) do Ch_
            for (Vl, D, ϵ, τ) in EXP_PARAMS
                put!(Ch_, (Vl, D, ϵ, τ))
            end
        end
        
        @threads for thid in 1:nths
            for (Vl, D, ϵ, τ) in Ch
                # LOAD
                status = MINDEX[:STATUS, Vl, D, ϵ, τ]
                status != :stst && continue

                # Dynamic
                M0 = idxdat([:M0], Vl, D, ϵ, τ)
                lock(WLOCK) do; isnothing(LP_cache) && (LP_cache = Dyn.vgvatp_cache(M0)) end
                z(vatp, vg) = LP_cache[vatp][vg][M0.obj_idx]
                vg(vatp, vg) = vg

                fX(vatp, vg) = M0.Xb[vatp][vg] / M0.X
                PX = Dyn.get_join(fX, M0)
                
                SPX = Dyn.entropy(PX)
                biom_avPX = Dyn.ave_over(z, PX)
                biom_vaPX = Dyn.ave_over((vatp_, vg_) -> (biom_avPX - vatp_)^2, PX)
                vg_avPX = Dyn.ave_over(vg, PX)
                vg_vaPX = Dyn.ave_over((vatp_, vg_) -> (vg_avPX - vg_)^2, PX) 

                lock(WLOCK) do
                    JDAT_[(:dyn, Vl, D, ϵ, τ, :S)] = SPX
                    JDAT_[(:dyn, Vl, D, ϵ, τ, :biom_av)] = biom_avPX
                    JDAT_[(:dyn, Vl, D, ϵ, τ, :biom_va)] = biom_vaPX
                    JDAT_[(:dyn, Vl, D, ϵ, τ, :vg_av)] = vg_avPX
                    JDAT_[(:dyn, Vl, D, ϵ, τ, :vg_va)] = vg_vaPX
                end
                
                # MaxEnt
                for MEmode in ME_MODELS
                    M = idxdat([MEmode, :M], Vl, D, ϵ, τ)

                    beta_biom = idxdat([MEmode, :beta_biom], Vl, D, ϵ, τ)
                    beta_vg = idxdat([MEmode, :beta_vg], Vl, D, ϵ, τ)
                    fME(vatp_, vg_) = exp((beta_biom * z(vatp_, vg_)) + (beta_vg * vg(vatp_, vg_)))
                    PME = Dyn.get_join(fME, M)
                    
                    SPME = Dyn.entropy(PME)
                    biom_avPME = Dyn.ave_over(z, PME)
                    biom_vaPME = Dyn.ave_over((vatp_, vg_) -> (biom_avPME - vatp_)^2, PME) 
                    vg_avPME = Dyn.ave_over(vg, PME)
                    vg_vaPME = Dyn.ave_over((vatp_, vg_) -> (vg_avPME - vg_)^2, PME) 

                    lock(WLOCK) do
                        JDAT_[(MEmode, Vl, D, ϵ, τ, :S)] = SPME
                        JDAT_[(MEmode, Vl, D, ϵ, τ, :biom_av)] = biom_avPME
                        JDAT_[(MEmode, Vl, D, ϵ, τ, :biom_va)] = biom_vaPME
                        JDAT_[(MEmode, Vl, D, ϵ, τ, :vg_av)] = vg_avPME
                        JDAT_[(MEmode, Vl, D, ϵ, τ, :vg_va)] = vg_vaPME

                        @info(string("Doing...", "-"^50), 
                            MEmode, (Vl, D, ϵ, τ), 
                            (SPX, SPME), 
                            (biom_avPX, biom_avPME), 
                            (biom_vaPX, biom_vaPME), 
                            (vg_avPX, vg_avPME),
                            (vg_vaPX, vg_vaPME)
                        ); println()
                    end
                end # for MEmode
            end # for Ch
        end # for thid

        return JDAT_
    end # load_cache do

    # ---------------------------------------------------------------
    # Plots
    let
        fontsize = 13
        sparams = (;
            # dpi = 1000,
            thickness_scaling = 1.3, 
            xguidefontsize = fontsize, yguidefontsize = fontsize
        )

        
        function datf(src, Vl, D, ϵ, τ, dkey)
            if dkey in [:S, :biom_av, :vg_av]
                return JDAT[(src, Vl, D, ϵ, τ, dkey)]
            elseif dkey == :vg_va
                std = sqrt(JDAT[(src, Vl, D, ϵ, τ, :vg_va)])
                # av = JDAT[(:dyn, Vl, D, ϵ, τ, :vg_va)]
                av = 1.0
                return std/av
            elseif dkey == :biom_va
                std = sqrt(JDAT[(src, Vl, D, ϵ, τ, :biom_va)])
                # av = JDAT[(:dyn, Vl, D, ϵ, τ, :biom_va)]
                av = 1.0
                return std/av
            end
        end

        for dkey in [
                :S, :biom_av, :biom_va, :vg_av, :vg_va
            ]
            ps_pool = Dict()
            vals_pool = Dict()
            for (Vl, D, ϵ, τ) in EXP_PARAMS

                # ϵ < 0.1 && continue

                # LOAD
                status = MINDEX[:STATUS, Vl, D, ϵ, τ]
                status != :stst && continue

                # valPX = f(JDAT[(:dyn, Vl, D, ϵ, τ, dkey)])
                valPX = datf(:dyn, Vl, D, ϵ, τ, dkey)

                for MEmode in ME_MODELS

                    p = get!(ps_pool, MEmode, 
                        plot(;
                            title = string(MEmode, " ($dkey)"),
                            xlabel = MODEL_LABELS[:dyn], ylabel = MODEL_LABELS[MEmode]
                        )
                    )

                    # valPME = f(JDAT[(MEmode, Vl, D, ϵ, τ, dkey)])
                    valPME = datf(MEmode, Vl, D, ϵ, τ, dkey)

                    color = ES_COLORS[ϵ]
                    marker = (8, MODEL_MARKERS[MEmode])
                    scatter!(p, [valPX], [valPME]; 
                        label = "", color, marker, alpha = 0.7, 
                    )

                    vals = get!(vals_pool, MEmode, [])
                    push!(vals, valPME, valPX)
                end
                
            end

            ps = Plots.Plot[]
            for (MEmode, p) in ps_pool
                vals = sort!(vals_pool[MEmode])
                m = (maximum(vals) - minimum(vals)) * 0.1
                plot!(p, vals, vals; label = "", 
                    lw = 3, alpha = 0.7, ls = :dash, 
                    xlim = [minimum(vals) - m, maximum(vals) + m],
                    ylim = [minimum(vals) - m, maximum(vals) + m],
                    sparams...
                )
                # mysavefig(p, string(dkey, "_corrs"); MEmode)
                push!(ps, p)
            end
            mysavefig(ps, string(dkey, "_corrs"))
        end
    end

end
  