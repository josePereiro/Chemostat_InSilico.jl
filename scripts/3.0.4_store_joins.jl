import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic
        
    using ProgressMeter
    using Plots.Measures

    using Plots
    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    using Base.Threads

    UJL.set_cache_dir(Dyn.cachedir())
end

# ----------------------------------------------------------------------------
# MINDEX
MINDEX_FILE = Dyn.procdir("marg_dat_index.bson")
MINDEX = UJL.load_data(MINDEX_FILE)
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...)
marg_dat(dk, indexks...) = Dyn.idxdat(MINDEX, dk, indexks...)

# ----------------------------------------------------------------------------
let
    
    ME_MODELS = [
        :ME_Z_OPEN_G_OPEN, 
        :ME_FULL_POLYTOPE,
        :ME_Z_EXPECTED_G_EXPECTED,
        :ME_Z_EXPECTED_G_BOUNDED,
        :ME_Z_FIXXED_G_BOUNDED, 
    ]
    
    # ---------------------------------------------------------------
    # COLLECT DATA
    # store collected data from join distributions
    nths = 3 # working threads    
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
            M0 = marg_dat([:M0], Vl, D, ϵ, τ)
            
            lock(WLOCK) do; isnothing(LP_cache) && (LP_cache = Dyn.vgvatp_cache(M0)) end
            z(vatp, vg) = LP_cache[vatp][vg][M0.obj_idx]
            vg(vatp, vg) = vg
            vatp(vatp, vg) = vatp

            fX(vatp, vg) = M0.Xb[vatp][vg] / M0.X
            PX = Dyn.get_join(fX, M0)
            
            SPX = Dyn.entropy(PX)
            biom_avPX = Dyn.ave_over(z, PX)
            biom_vaPX = Dyn.ave_over((vatp_, vg_) -> (biom_avPX - z(vatp_, vg_))^2, PX)
            vatp_avPX = Dyn.ave_over(vatp, PX)
            vatp_vaPX = Dyn.ave_over((vatp_, vg_) -> (vatp_avPX - vatp_)^2, PX) 
            vg_avPX = Dyn.ave_over(vg, PX)
            vg_vaPX = Dyn.ave_over((vatp_, vg_) -> (vg_avPX - vg_)^2, PX) 

            let
                JDAT = Dict()
                JDAT[:P] = PX
                JDAT[:S] = SPX
                JDAT[:biom_av] = biom_avPX
                JDAT[:biom_va] = biom_vaPX
                JDAT[:vatp_av] = vatp_avPX
                JDAT[:vatp_va] = vatp_vaPX
                JDAT[:vg_av] = vg_avPX
                JDAT[:vg_va] = vg_vaPX

                mode = :dyn
                @info(string("Done  ", "-"^50), 
                    mode, (Vl, D, ϵ, τ)
                ); println()

                JDAT_CID = (mode, Vl, D, ϵ, τ)
                UJL.exist_cache(JDAT_CID) || UJL.save_cache(JDAT_CID, JDAT; verbose = false)
            end

            # MaxEnt
            for MEmode in ME_MODELS

                M = marg_dat([MEmode, :M], Vl, D, ϵ, τ)

                beta_biom = marg_dat([MEmode, :beta_biom], Vl, D, ϵ, τ)
                beta_vg = marg_dat([MEmode, :beta_vg], Vl, D, ϵ, τ)
                fME(vatp_, vg_) = exp((beta_biom * z(vatp_, vg_)) + (beta_vg * vg(vatp_, vg_)))
                PME = Dyn.get_join(fME, M)
                
                SPME = Dyn.entropy(PME)
                biom_avPME = Dyn.ave_over(z, PME)
                biom_vaPME = Dyn.ave_over((vatp_, vg_) -> (biom_avPME - z(vatp_, vg_))^2, PME) 
                vatp_avPME = Dyn.ave_over(vatp, PME)
                vatp_vaPME = Dyn.ave_over((vatp_, vg_) -> (vatp_avPME - vatp_)^2, PME) 
                vg_avPME = Dyn.ave_over(vg, PME)
                vg_vaPME = Dyn.ave_over((vatp_, vg_) -> (vg_avPME - vg_)^2, PME) 

                let
                    JDAT = Dict()
                    JDAT[:P] = PME
                    JDAT[:S] = SPME
                    JDAT[:biom_av] = biom_avPME
                    JDAT[:biom_va] = biom_vaPME
                    JDAT[:vatp_av] = vatp_avPME
                    JDAT[:vatp_va] = vatp_vaPME
                    JDAT[:vg_av] = vg_avPME
                    JDAT[:vg_va] = vg_vaPME

                    @info(string("Done...", "-"^50), 
                        MEmode, (Vl, D, ϵ, τ), 
                        (SPX, SPME), 
                        (biom_avPX, biom_avPME), 
                        (biom_vaPX, biom_vaPME), 
                        (vg_avPX, vg_avPME),
                        (vg_vaPX, vg_vaPME)
                    ); println()

                    JDAT_CID = (MEmode, Vl, D, ϵ, τ)
                    UJL.exist_cache(JDAT_CID) || UJL.save_cache(JDAT_CID, JDAT; verbose = false)
                end
            end # for MEmode

        end # for Ch
    end # for thid
end