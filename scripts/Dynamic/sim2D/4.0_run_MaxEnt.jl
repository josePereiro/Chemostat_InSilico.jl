@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        SimD2, run_simD2!, check_stst, hist, 
        n_ave_conv!, tail_ave,
        plotsdir, lglob, sglob, Container, vec!, Vi, 
        normalizeP!, 
        save_simdat, load_simdat, simdat_file,
        set_status, get_status,
        load_batch, save_batch,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS, collect_ts,
        STST, is_steady
        
    import UtilsJL
    const PltU = UtilsJL.PlotsUtils
    const Ass = UtilsJL.ProjAssistant
    const GU = UtilsJL.GeneralUtils

    using Plots
    import GR
    !isinteractive() && GR.inline("png")

    using Base.Threads
    using ProgressMeter
    using BenchmarkTools
    using ExtractMacro

    import DataFileNames
    const DFN = DataFileNames

    Ass.set_verbose(false)
end

## ------------------------------------------------------
# globlas
# batch = (; push_frec, Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts)
Ds, ϵs, cgs, simid = Dyn.lglob(:Ds, :ϵs, :cgs, :SimD2Id)
iter = collect(Iterators.product(Ds, ϵs, cgs))
simparams = filter(iter) do tp
    D, ϵ, cg = tp
    isinf(cg) && return false
    status = Dyn.get_status(simid, (;D, ϵ, cg))
    (status != STST_SIM_STATUS) && return false
    return true
end;

# ME_AVZ_EXPECTED_AVUG_BOUNDED
# ME_AVZ_EXPECTED_AVUG_EXPECTED

## ----------------------------------------------------------------------------
let

    WLOCK = ReentrantLock()
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    sim_tot = length(simparams)
    # @threads 
    for (simi, (D, ϵ, cg)) in collect(enumerate(simparams))

        # Test
        D=3.86e-01 
        cg=1.50e+01 
        ϵ=1.00e-02

        thid = threadid()

        simparams = (;D, ϵ, cg)

        # ---------------------------------------------------------------
        # Sim globals
        S = Dyn.load_simdat(simid, simparams)
        X = S.X 
        
        # ---------------------------------------------------------------
        # Space globals
        V = S.V
        Vz = Vi(V, :z)
        z0 = sum(Vz) / length(Vz)
        Vug = Vi(V, :ug)
        ug0 = sum(Vug) / length(Vug)
        z_avPX = S.z_av
        ug_avPX = S.ug_av
        cgD_X = S.cg * S.D / S.X
        S = nothing; GC.gc()

        # ---------------------------------------------------------------
        # MaxEnt globals
        PME = ones(length(V))
        z_avPME = 0.0
        ug_avPME = 0.0
        normalizeP!(PME)
        z_beta = 0.0
        ug_beta = 0.0
        ug_avPME_at_b0 = 0.0
        z_betas = Float64[]
        ug_betas = Float64[]


        ## -----------------------------------------------------------
        topiter = 1
        top_maxiter = 100
        conv = false

        ## -----------------------------------------------------------
        # GD globals
        gd_maxiter = 1000
        gdth = 1e-2
        gd_err = 0.0
        verb_frec = 100
        gditer = 0
        
        ## -----------------------------------------------------------
        # beta stst
        stst_th, stst_w = 1e-2, 8
        nbuffs = length((:z_beta, :ug_beta))
        beta_stst = STST(nbuffs, stst_w)
        
        ## -----------------------------------------------------------
        # Info funtions
        function gd_ug_info() 
            lock(WLOCK) do
                @info("ug grad descent ",
                    (simi, sim_tot),
                    topiter, gditer,  
                    MEmode, 
                    (ug_avPME, cgD_X), 
                    gd_err, z_beta, thid
                ); println()
            end
        end

        function gd_z_info()
            lock(WLOCK) do
                @info("z grad Descent ",
                    (simi, sim_tot),
                    topiter, gditer, 
                    MEmode, 
                    (z_avPME, z_avPX), 
                    gd_err, z_beta, thid
                ); println()
            end
        end

        function finished_round_info()
            lock(WLOCK) do
                @info("Finished Round",
                    (simi, sim_tot),
                    topiter, conv,
                    MEmode,
                    (z_beta, ug_beta),
                    (ug_avPME, cgD_X),
                    (ug_avPME_at_b0, cgD_X),
                    (z_avPME, z_avPX),
                    thid
                ); println()
            end
        end

        ## -----------------------------------------------------------
        while true

            ## -----------------------------------------------------------
            # find z_beta
            # Gradient descent
            target = z_avPX
            x0 = z_beta
            maxΔx = max(100.0, abs(z_beta) * 0.1)
            x1 = x0 + maxΔx * 0.1

            function gd_a_upfun!(gdmodel)

                gditer = gdmodel.iter
                z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

                PME .= exp.(z_beta .* (Vz .- z0)) .* exp.(ug_beta .* (Vug .- ug0))
                normalizeP!(PME)

                z_avPME = sum(PME .* Vz)
                
                gd_err = gdmodel.ϵi

                return z_avPME
            end

            gdmodel = UtilsJL.SimulationUtils.grad_desc(gd_a_upfun!; gdth, 
                target, x0, x1, maxΔx, 
                maxiter = gd_maxiter, verbose = false
            )
            z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

            gd_z_info()
            
            ## -----------------------------------------------------------
            # Check balance at ug_beta = 0.0
            ug_beta0 = 0.0
            PME .= exp.(z_beta .* (Vz .- z0)) .* exp.(ug_beta0 .* (Vug .- ug0))
            normalizeP!(PME)
            ug_avPME_at_b0 = sum(PME .* Vug)
            ug_valid_at_b0 = (ug_avPME_at_b0 <= cgD_X)

            ## -----------------------------------------------------------
            # if not valid, move ug_beta
            if ug_valid_at_b0
                ug_beta = 0.0
            else
                # Gradient descent
                target = cgD_X * (1.0 - gdth)
                x0 = ug_beta
                maxΔx = max(10.0, abs(ug_beta) * 0.1)
                x1 = x0 + maxΔx * 0.1
                
                function ug_ug_upfun!(gdmodel)

                    gditer = gdmodel.iter
                    ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

                    PME .= exp.(z_beta .* (Vz .- z0)) .* exp.(ug_beta .* (Vug .- ug0))
                    normalizeP!(PME)
                    ug_avPME = sum(PME .* Vug)
                    
                    gd_err = gdmodel.ϵi

                    return ug_avPME 
                end

                gdmodel = UtilsJL.SimulationUtils.grad_desc(ug_ug_upfun!; gdth, 
                    target, x0, x1, maxΔx, 
                    maxiter = gd_maxiter, verbose = false
                )
                ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
            end

            gd_ug_info()

            push!(z_betas, z_beta); push!(ug_betas, ug_beta)
            push!(beta_stst, z_beta, ug_beta)

            conv = ug_valid_at_b0 || is_steady(beta_stst, stst_th)
            
            finished_round_info()

            conv && break
            
            topiter += 1
            topiter > top_maxiter && break
        end

        # save
        dat = (;PME, z_beta, ug_beta, z_avPME, ug_avPME)
        Dyn.save_maxent(dat, simid, MEmode, simparams)

        break
    end # for (simi, (D, ϵ, cg)) 

end

## ----------------------------------------------------------------------------
# function run_ME!(S, D, X, MEmode)
    
#     thid = threadid()

#     ## -----------------------------------------------------------
#     # Globals
#     cgD_X = M.cg * M.D/ M.X
#     biom_idx = M.obj_idx
#     vl_idx = M.vl_idx
#     ug_idx = M.ug_idx
    
#     ## -----------------------------------------------------------
#     # G BOUNDING

#     ## -----------------------------------------------------------
#     is_bounded = MEmode in [ME_Z_OPEN_G_BOUNDED, ME_Z_EXPECTED_G_BOUNDED, ME_Z_FIXXED_G_BOUNDED]
#     is_bounded && let
#         # Fix av_ug
#         net = M.net
#         net.ub[ug_idx] = min(M.Vg, cgD_X)
#         net.ub[vl_idx] = min(M.Vl, cgD_X)
#         L, U = Dyn.fva(net)
#         net.lb .= L; net.ub .= U
#     end

#     ## -----------------------------------------------------------
#     # Z FIXXED

#     ## -----------------------------------------------------------
#     is_fixxed = MEmode in [ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED]
#     is_fixxed && let
#         # Fix biomass to observable
#         net = M.net
#         net.ub[biom_idx] = z_avPX * (1.0 + δμ)
#         net.lb[biom_idx] = z_avPX * (1.0 - δμ)
#         L, U = Dyn.fva(net)
#         net.lb .= L; net.ub .= U
#     end

#     ## -----------------------------------------------------------
#     # EXPECTED
    
#     ## -----------------------------------------------------------
#     # Globals
#     PME = nothing
#     z_beta = 0.0
#     ug_beta = 0.0
#     maxiter = 800
#     verb_frec = 50.0
#     gdth = 1e-2

#     ## -----------------------------------------------------------
#     # Z EXPECTED

#     ## -----------------------------------------------------------
#     is_zexpected = MEmode in [ME_Z_EXPECTED_G_OPEN, ME_Z_EXPECTED_G_BOUNDED]
#     is_zexpected && let
#         # Gradient descent
#         target = z_avPX
#         x0 = 1.5e2
#         x1 = x0 * 0.9
#         maxΔx = 100.0
#         gdit = 1

#         # grad desc
#         PME = Dyn.get_join(M)
#         function gd_a_upfun!(gdmodel)
            
#             z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

#             PME = Dyn.get_join!(M, PME) do vatp_, ug_
#                 exp(z_beta * z(vatp_, ug_))
#             end
#             z_avPME = Dyn.ave_over(z, PME)
            
#             biom_err = gdmodel.ϵi
#             show_info = gdit == 1 || rem(gdit, verb_frec) == 0 || 
#                 gdit == maxiter || biom_err < gdth
#             show_info && lock(WLOCK) do
#                 @info("Grad Descent ", 
#                     gdit, MEmode, 
#                     (z_avPX, z_avPME), 
#                     biom_err, z_beta, thid
#                 ); println()
#             end

#             gdit += 1
#             return z_avPME
#         end

#         gdmodel = UtilsJL.SimulationUtils.grad_desc(gd_a_upfun!; gdth, 
#             target, x0, x1, maxΔx, maxiter, 
#             verbose = false
#         )
#         z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
#     end

#     ## -----------------------------------------------------------
#     # Z AND G EXPECTED

#     ## -----------------------------------------------------------
#     is_zgexpected = MEmode == ME_Z_EXPECTED_G_EXPECTED
#     is_zgexpected && let
#         target = [z_avPX, ug_avPX]
#         x0 = [100.0, -10.0]
#         x1 = [101.0, -11.0]
#         maxΔx = [80.0, 30.0]
#         gdit = 1

#         PME = Dyn.get_join(M)
#         function gd_biom_vg!(gdmodel)

#             z_beta, ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
#             PME = Dyn.get_join!(M, PME) do vatp_, ug_
#                 exp(z_beta * z(vatp_, ug_) + ug_beta * vg(vatp_, ug_))
#             end
#             z_avPME = Dyn.ave_over(z, PME)
#             ug_avPME = Dyn.ave_over(vg, PME)
            
#             biom_err = abs(z_avPX - z_avPME)/z_avPX
#             ug_err = abs(ug_avPX - ug_avPME)/ug_avPX
#             gd_err = max(biom_err, ug_err)

#             show_info = gdit == 1 || rem(gdit, verb_frec) == 0 || 
#                 gdit == maxiter || gd_err < gdth
#             show_info && lock(WLOCK) do
#                 @info("Grad Descent ", 
#                     gdit, MEmode, 
#                     (z_avPX, z_avPME),
#                     (ug_avPX, ug_avPME),
#                     gd_err, z_beta, ug_beta, 
#                     thid
#                 ); println()
#             end

#             gdit += 1
#             return [z_avPME, ug_avPME]
#         end

#         gdmodel = UJL.grad_desc_vec(gd_biom_vg!; 
#             target, x0, x1, maxΔx, gdth, maxiter, 
#             verbose = false
#         )
#         z_beta, ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
#     end

#     ## -----------------------------------------------------------
#     is_full = MEmode == ME_FULL_POLYTOPE
#     is_full && let
#         topiter = 1
#         PME = Dyn.get_join(M)
#         stth, stw = 0.1, 8
#         z_betas, ug_betas = [z_beta], [ug_beta]
#         z_avPME, ug_avPME = 0.0, 0.0

#         while true

#             ## -----------------------------------------------------------
#             # find z_beta
#             # Gradient descent
#             target = z_avPX
#             x0 = z_beta
#             maxΔx = max(100.0, abs(z_beta) * 0.1)
#             x1 = x0 + maxΔx * 0.1

#             function gd_a_upfun!(gdmodel)

#                 gditer = gdmodel.iter
#                 z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

#                 PME = Dyn.get_join!(M, PME) do vatp_, ug_
#                     exp(z_beta * z(vatp_, ug_) + ug_beta * vg(vatp_, ug_))
#                 end
#                 z_avPME = Dyn.ave_over(z, PME)
                
#                 gd_err = gdmodel.ϵi
#                 show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
#                     gditer == gd_maxiter || gd_err < gdth
#                 show_info && begin
#                     @info("z grad Descent ", 
#                         topiter, gditer, 
#                         MEmode, 
#                         (z_avPX, z_avPME), 
#                         gd_err, z_beta, thid
#                     ); println()
#                 end

#                 return z_avPME 
#             end

#             gdmodel = UtilsJL.SimulationUtils.grad_desc(gd_a_upfun!; gdth, 
#                 target, x0, x1, maxΔx, 
#                 maxiter, verbose = false
#             )
#             z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
            
#             ## -----------------------------------------------------------
#             # Check balance at ug_beta = 0.0
#             PME = Dyn.get_join!(M, PME) do vatp_, ug_
#                 ug_beta_ = 0.0
#                 exp(z_beta * z(vatp_, ug_) + ug_beta_ * vg(vatp_, ug_))
#             end
#             ug_avPME_at_b0 = Dyn.ave_over(vg, PME)
#             ug_valid_at_b0 = ug_avPME_at_b0 <= cgD_X

#             ## -----------------------------------------------------------
#             # if not valid, move ug_beta
#             if ug_valid_at_b0
#                 ug_beta = 0.0
#             else
#                 # Gradient descent
#                 target = cgD_X * (1.0 - gdth)
#                 x0 = ug_beta
#                 maxΔx = max(10.0, abs(ug_beta) * 0.1)
#                 x1 = x0 + maxΔx * 0.1
                
#                 function ug_ug_upfun!(gdmodel)

#                     gditer = gdmodel.iter
#                     ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

#                     PME = Dyn.get_join!(M, PME) do vatp_, ug_
#                         exp(z_beta * z(vatp_, ug_) + ug_beta * vg(vatp_, ug_))
#                     end
#                     ug_avPME = Dyn.ave_over(vg, PME)
                    
#                     gd_err = gdmodel.ϵi
#                     show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
#                         gditer == gd_maxiter || gd_err < gdth
#                     show_info && begin
#                         @info("vg grad descent ", 
#                             topiter, gditer,  
#                             MEmode, 
#                             (cgD_X, ug_avPME), 
#                             gd_err, z_beta, thid
#                         ); println()
#                     end

#                     return ug_avPME 
#                 end

#                 gdmodel = UtilsJL.SimulationUtils.grad_desc(ug_ug_upfun!; gdth, 
#                     target, x0, x1, maxΔx, 
#                     maxiter, verbose = false
#                 )
#                 ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
#             end

#             push!(z_betas, z_beta); push!(ug_betas, ug_beta)

#             conv = (ug_valid_at_b0 && gdmodel.ϵi < gdth) || 
#                 (UJL.is_stationary(z_betas, stth, stw) && UJL.is_stationary(ug_betas, stth, stw))
            
#             @info("Finished Round", 
#                 topiter, conv,
#                 MEmode,
#                 (z_beta, ug_beta),
#                 (ug_avPME, cgD_X),
#                 (ug_avPME_at_b0, cgD_X),
#                 (z_avPME, z_avPX),
#                 thid
#             ); println()

#             conv && break
            
#             topiter += 1
#             topiter > maxiter && break
#         end
#     end

#     ## -----------------------------------------------------------
#     is_moving = MEmode == ME_Z_EXPECTED_G_MOVING
#     is_moving && let

#         # init globals
#         Δstep = 0.5
#         maxiter = 500
#         z_avPME = 0.0
#         biom_err = Inf
        
#         ## -----------------------------------------------------------
#         for rit in 1:maxiter

#             ## -----------------------------------------------------------
#             # Find biom beta
#             target = z_avPX
#             x0 = z_beta
#             x1 = x0 * 0.9
#             maxΔx = 100.0
        
#             # grad desc
#             it = 1
#             PME = Dyn.get_join(M)
#             function gd_a_upfun!(gdmodel)
#                 beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

#                 PME = Dyn.get_join!(M, PME) do vatp_, ug_
#                     exp(beta * z(vatp_, ug_))
#                 end
#                 z_avPME = Dyn.ave_over(z, PME)
                
#                 biom_err = gdmodel.ϵi
#                 show_info = it == 1 || rem(it, verb_frec) == 0 || 
#                     it == maxiter || biom_err < gdth
#                 show_info && lock(WLOCK) do
#                     @info("Grad Descent ", 
#                         it, MEmode, 
#                         (z_avPX, z_avPME), 
#                         biom_err, beta, thid
#                     ); println()
#                 end

#                 it += 1
#                 return z_avPME
#             end

#             gdmodel = UtilsJL.SimulationUtils.grad_desc(gd_a_upfun!; gdth, 
#                 target, x0, x1, maxΔx, maxiter, 
#                 verbose = false
#             )
#             z_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

#             ## -----------------------------------------------------------
#             # move vg
#             ug_avPME = Dyn.ave_over(vg, PME)
        
#             Δ = cgD_X - ug_avPME
#             net = M.net
#             net.ub[ug_idx] = max(net.lb[ug_idx], 
#                 min(M.Vg, net.ub[ug_idx] + Δstep * Δ)
#             )
#             L, U = Dyn.fva(net)
#             net.lb .= L; net.ub .= U
#             ug_ub = net.ub[ug_idx]
        
#             valid_ug_avPME = Δ >= 0

#             ## -----------------------------------------------------------
#             lock(WLOCK) do 
#                 @info("End round", 
#                     rit, z_beta, 
#                     biom_err, valid_ug_avPME, 
#                     (M.Vg, ug_ub),
#                     (cgD_X, ug_avPME),
#                     thid
#                 ); println()
#             end

#             conv = valid_ug_avPME && biom_err < gdth
#             conv && break

#         end # for rit in 1:maxiter
#     end

#     ## -----------------------------------------------------------
#     # DO ME AND MARGINALIZE
#     MEMs = Dyn.get_marginals(M; δ, LP_cache, verbose = false) do vatp_, ug_
#         exp((z_beta * z(vatp_, ug_)) + (ug_beta * vg(vatp_, ug_)))
#     end
#     return PME, MEMs, z_beta, ug_beta
# end
