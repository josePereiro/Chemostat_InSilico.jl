## ----------------------------------------------------------------------------
function maxent(V, D, cgD_X; 
        # Top loop
        top_maxiter = 100,

        # GD globals
        gd_maxiter = 1000,
        gdth = 1e-2,
        verbose = false,

        # stst
        stst_th = 1e-2,
        stst_w = 10,

        check_b0 = true

    )

    # ---------------------------------------------------------------
    # Space globals
    Vz = Vi(V, :z)
    z0 = sum(Vz) / length(Vz)
    Vz0 = Vz .- z0
    Vug = Vi(V, :ug)
    ug0 = sum(Vug) / length(Vug)
    Vug0 = Vug .- ug0

    # ---------------------------------------------------------------
    # Alg globals
    PME = ones(length(V))
    normalizeP!(PME)
    Paux = similar(PME)
    I = eachindex(PME)
    z_avPME = 0.0
    ug_avPME = 0.0
    z_beta = 0.0
    ug_beta = 0.0
    ug_avPME_at_b0 = 0.0
    ug_valid_at_b0 = false
    z_betas = Float64[]
    ug_betas = Float64[]

    ## -----------------------------------------------------------
    # top loop globals
    topiter = 1
    conv = false
    thid = threadid()

    ## -----------------------------------------------------------
    # GD globals
    gd_err = 0.0
    gd_iter = 0

    ## -----------------------------------------------------------
    # beta stst
    nbuffs = length((:z_beta, :ug_beta))
    beta_stst = STST(nbuffs, stst_w)

    ## -----------------------------------------------------------
    # Info funtions
    function gd_ug_info() 
        @info("ug grad descent ",
            topiter, gd_iter,   
            (ug_avPME, cgD_X), 
            gd_err, ug_beta, thid
        )
    end

    function gd_z_info()
        @info("z grad Descent ",
            topiter, gd_iter, 
            (z_avPME, D), 
            gd_err, z_beta, thid
        )
    end

    function finished_round_info()
        @info("Finished Round",
            topiter, conv,
            (z_beta, ug_beta),
            (ug_avPME, cgD_X),
            (ug_avPME_at_b0, cgD_X),
            (z_avPME, D),
            thid
        ); 
    end

    ## -----------------------------------------------------------
    while true

        ## -----------------------------------------------------------
        # find z_beta
        # Gradient descent
        target = D
        x0 = z_beta
        maxΔx = max(100.0, abs(z_beta) * 0.1)
        x1 = x0 + maxΔx * 0.1
        gd_iter = 1

        function gd_z_upfun!(gdmodel)

            z_beta = SimTools.gd_value(gdmodel)

            # aux
            (gd_iter == 1) && (Paux .= exp.(ug_beta .* Vug0))

            PME .= exp.(z_beta .* Vz0) .* Paux
            normalizeP!(PME)

            z_avPME = 0.0
            @inbounds for i in I
                z_avPME += PME[i] * Vz[i]
            end
            
            gd_iter += 1
            gd_err = gdmodel.ϵi

            return z_avPME
        end

        gdmodel = SimTools.grad_desc(gd_z_upfun!; gdth, 
            target, x0, x1, maxΔx, 
            maxiter = gd_maxiter, verbose = false
        )
        z_beta = SimTools.gd_value(gdmodel)

        verbose && gd_z_info()
        
        ## -----------------------------------------------------------
        # TODO: prove all this
        # Check balance at ug_beta = 0.0 the first time
        if check_b0 && (topiter == 1)
            ug_beta0 = 0.0
            PME .= exp.(z_beta .* Vz0) .* exp.(ug_beta0 .* Vug0)
            normalizeP!(PME)
            ug_avPME_at_b0 = sum(PME .* Vug)
            ug_valid_at_b0 = (abs(ug_avPME_at_b0) <= abs(cgD_X))
        end

        ## -----------------------------------------------------------
        # if not valid, move ug_beta
        if ug_valid_at_b0
            ug_beta = 0.0
            ug_avPME = ug_avPME_at_b0
        else
            # Gradient descent
            target = cgD_X * (1.0 - gdth)
            x0 = ug_beta
            maxΔx = max(10.0, abs(ug_beta) * 0.1)
            x1 = x0 + maxΔx * 0.1
            gd_iter = 1
            
            function ug_ug_upfun!(gdmodel)

                ug_beta = SimTools.gd_value(gdmodel)

                (gd_iter == 1) && (Paux .= exp.(z_beta .* Vz0))

                PME .= Paux .* exp.(ug_beta .* Vug0)
                normalizeP!(PME)

                ug_avPME = 0.0
                @inbounds for i in I
                    ug_avPME += PME[i] * Vug[i]
                end
                
                gd_iter += 1
                gd_err = gdmodel.ϵi

                return ug_avPME 
            end

            gdmodel = SimTools.grad_desc(ug_ug_upfun!; gdth, 
                target, x0, x1, maxΔx, 
                maxiter = gd_maxiter, verbose = false
            )
            ug_beta = SimTools.gd_value(gdmodel)
        end

        verbose && gd_ug_info()

        push!(z_betas, z_beta); push!(ug_betas, ug_beta)
        push!(beta_stst, z_beta, ug_beta)

        conv = ug_valid_at_b0 || is_steady(beta_stst, stst_th)
        
        verbose && finished_round_info()

        conv && break
        
        topiter += 1
        topiter > top_maxiter && break
    end # while true

    status = conv ? :conv : :maxiter
    return (;PME, z_beta, ug_beta, z_avPME, ug_avPME, status)
    
end

match_moms(V, D, cgD_X; kwargs...) = maxent(V, D, cgD_X; kwargs..., check_b0 = false)