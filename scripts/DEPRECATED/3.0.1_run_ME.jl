## ----------------------------------------------------------------------------
function run_ME!(M, MEmode; LP_cache, δ, δμ, biom_avPX, vg_avPX)
    
    thid = threadid()

    ## -----------------------------------------------------------
    # Globals
    cgD_X = M.cg * M.D/ M.X
    biom_idx = M.obj_idx
    vl_idx = M.vl_idx
    vg_idx = M.vg_idx
    z(vatp, vg) = LP_cache[vatp][vg][biom_idx]
    vg(vatp, vg) = vg
    
    ## -----------------------------------------------------------
    # G BOUNDING

    ## -----------------------------------------------------------
    is_bounded = MEmode in [ME_Z_OPEN_G_BOUNDED, ME_Z_EXPECTED_G_BOUNDED, ME_Z_FIXXED_G_BOUNDED]
    is_bounded && let
        # Fix av_ug
        net = M.net
        net.ub[vg_idx] = min(M.Vg, cgD_X)
        net.ub[vl_idx] = min(M.Vl, cgD_X)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    ## -----------------------------------------------------------
    # Z FIXXED

    ## -----------------------------------------------------------
    is_fixxed = MEmode in [ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED]
    is_fixxed && let
        # Fix biomass to observable
        net = M.net
        net.ub[biom_idx] = biom_avPX * (1.0 + δμ)
        net.lb[biom_idx] = biom_avPX * (1.0 - δμ)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    ## -----------------------------------------------------------
    # EXPECTED
    
    ## -----------------------------------------------------------
    # Globals
    PME = nothing
    beta_biom = 0.0
    beta_vg = 0.0
    maxiter = 800
    verb_frec = 50.0
    gdth = 1e-2

    ## -----------------------------------------------------------
    # Z EXPECTED

    ## -----------------------------------------------------------
    is_zexpected = MEmode in [ME_Z_EXPECTED_G_OPEN, ME_Z_EXPECTED_G_BOUNDED]
    is_zexpected && let
        # Gradient descent
        target = biom_avPX
        x0 = 1.5e2
        x1 = x0 * 0.9
        maxΔx = 100.0
        gdit = 1

        # grad desc
        PME = Dyn.get_join(M)
        function gd_biom!(gdmodel)
            
            beta_biom = UJL.gd_value(gdmodel)

            PME = Dyn.get_join!(M, PME) do vatp_, vg_
                exp(beta_biom * z(vatp_, vg_))
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

        gdmodel = UJL.grad_desc(gd_biom!; gdth, 
            target, x0, x1, maxΔx, maxiter, 
            verbose = false
        )
        beta_biom = UJL.gd_value(gdmodel)
    end

    ## -----------------------------------------------------------
    # Z AND G EXPECTED

    ## -----------------------------------------------------------
    is_zgexpected = MEmode == ME_Z_EXPECTED_G_EXPECTED
    is_zgexpected && let
        target = [biom_avPX, vg_avPX]
        x0 = [100.0, -10.0]
        x1 = [101.0, -11.0]
        maxΔx = [80.0, 30.0]
        gdit = 1

        PME = Dyn.get_join(M)
        function gd_biom_vg!(gdmodel)

            beta_biom, beta_vg = UJL.gd_value(gdmodel)
            PME = Dyn.get_join!(M, PME) do vatp_, vg_
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

        gdmodel = UJL.grad_desc_vec(gd_biom_vg!; 
            target, x0, x1, maxΔx, gdth, maxiter, 
            verbose = false
        )
        beta_biom, beta_vg = UJL.gd_value(gdmodel)
    end

    ## -----------------------------------------------------------
    is_full = MEmode == ME_FULL_POLYTOPE
    is_full && let
        topiter = 1
        PME = Dyn.get_join(M)
        stth, stw = 0.1, 8
        betas_biom, betas_vg = [beta_biom], [beta_vg]
        biom_avPME, vg_avPME = 0.0, 0.0

        while true

            ## -----------------------------------------------------------
            # find z_beta
            # Gradient descent
            target = biom_avPX
            x0 = beta_biom
            maxΔx = max(100.0, abs(beta_biom) * 0.1)
            x1 = x0 + maxΔx * 0.1

            function gd_biom!(gdmodel)

                gditer = gdmodel.iter
                beta_biom = UJL.gd_value(gdmodel)

                PME = Dyn.get_join!(M, PME) do vatp_, vg_
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

            gdmodel = UJL.grad_desc(gd_biom!; gdth, 
                target, x0, x1, maxΔx, 
                maxiter, verbose = false
            )
            beta_biom = UJL.gd_value(gdmodel)
            
            ## -----------------------------------------------------------
            # Check balance at beta_vg = 0.0
            PME = Dyn.get_join!(M, PME) do vatp_, vg_
                beta_vg_ = 0.0
                exp(beta_biom * z(vatp_, vg_) + beta_vg_ * vg(vatp_, vg_))
            end
            vg_avPME_at_b0 = Dyn.ave_over(vg, PME)
            vg_valid_at_b0 = vg_avPME_at_b0 <= cgD_X

            ## -----------------------------------------------------------
            # if not valid, move beta_vg
            if vg_valid_at_b0
                beta_vg = 0.0
            else
                # Gradient descent
                target = cgD_X * (1.0 - gdth)
                x0 = beta_vg
                maxΔx = max(10.0, abs(beta_vg) * 0.1)
                x1 = x0 + maxΔx * 0.1
                
                function gd_vg!(gdmodel)

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

                gdmodel = UJL.grad_desc(gd_vg!; gdth, 
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
    is_moving = MEmode == ME_Z_EXPECTED_G_MOVING
    is_moving && let

        # init globals
        Δstep = 0.5
        maxiter = 500
        biom_avPME = 0.0
        biom_err = Inf
        
        ## -----------------------------------------------------------
        for rit in 1:maxiter

            ## -----------------------------------------------------------
            # Find biom beta
            target = biom_avPX
            x0 = beta_biom
            x1 = x0 * 0.9
            maxΔx = 100.0
        
            # grad desc
            it = 1
            PME = Dyn.get_join(M)
            function gd_biom!(gdmodel)
                beta = UJL.gd_value(gdmodel)

                PME = Dyn.get_join!(M, PME) do vatp_, vg_
                    exp(beta * z(vatp_, vg_))
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

            gdmodel = UJL.grad_desc(gd_biom!; gdth, 
                target, x0, x1, maxΔx, maxiter, 
                verbose = false
            )
            beta_biom = UJL.gd_value(gdmodel)

            ## -----------------------------------------------------------
            # move vg
            vg_avPME = Dyn.ave_over(vg, PME)
        
            Δ = cgD_X - vg_avPME
            net = M.net
            net.ub[vg_idx] = max(net.lb[vg_idx], 
                min(M.Vg, net.ub[vg_idx] + Δstep * Δ)
            )
            L, U = Dyn.fva(net)
            net.lb .= L; net.ub .= U
            vg_ub = net.ub[vg_idx]
        
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
            conv && break

        end # for rit in 1:maxiter
    end

    ## -----------------------------------------------------------
    # DO ME AND MARGINALIZE
    MEMs = Dyn.get_marginals(M; δ, LP_cache, verbose = false) do vatp_, vg_
        exp((beta_biom * z(vatp_, vg_)) + (beta_vg * vg(vatp_, vg_)))
    end
    return PME, MEMs, beta_biom, beta_vg
end
