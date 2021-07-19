## ------------------------------------------------------
let

    WLOCK = ReentrantLock()
    MEmode = "ME_AVZ_EXPECTED_AVUG_EXPECTED"

    sim_tot = length(simparams)
    # @threads 
    for (simi, (D, ϵ, cg)) in collect(enumerate(simparams))

        # Test
        D=2.71e-01 
        cg=1.50e+01 
        ϵ=1.51e-01

        thid = threadid()

        simparams = (;D, ϵ, cg)

        # Sim
        S = Dyn.load_simdat(simid, simparams)
        X = S.X 
        
        # Space
        V = S.V
        Vz = Vi(V, :z)
        z0 = sum(Vz) / length(Vz)
        Vug = Vi(V, :ug)
        ug0 = sum(Vug) / length(Vug)
        z_avPX = S.z_av
        ug_avPX = S.ug_av
        S = nothing; GC.gc()

        # ---------------------------------------------------------------
        # MaxEnt 
        # TODO: This is not converging
        
        PME = ones(length(V))
        z_avPME = 0.0
        ug_avPME = 0.0
        normalizeP!(PME)
        z_beta = 0.0
        ug_beta = 0.0

        # gd globals
        target = [z_avPX, ug_avPX]
        x1 = [z_beta, ug_beta]
        maxΔx = [10.0, 10.0]
        minΔx = [0.05, 0.05]
        x0 = x1 .+ maxΔx .* 0.01
        gdit = 1
        maxsimparams = 1_000
        verb_frec = 100
        gdth = 1e-3
        gd_err = 0.0
        
        # gd fun
        function gd_up_fun!(gdmodel) 

            z_beta, ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)
            PME .= exp.(z_beta .* (Vz .- z0)) .* exp.(ug_beta .* (Vug .- ug0))
            normalizeP!(PME)

            z_avPME = sum(PME .* Vz)
            ug_avPME = sum(PME .* Vug)

            gd_err = maximum(gdmodel.ϵi)

            show_info = gdit == 1
            show_info |= rem(gdit, verb_frec) == 0 
            show_info |= gdit == maxsimparams 
            show_info |= gd_err < gdth
            show_info = false
            show_info && lock(WLOCK) do
                @info("Grad Descent ",
                    gdit, MEmode,
                    (simi, sim_tot), 
                    (z_avPX, z_avPME),
                    (ug_avPX, ug_avPME),
                    gd_err, z_beta, ug_beta,
                    thid
                ); println()
            end
            gdit += 1

            ret = [z_avPME, ug_avPME]
            any(isnan.(ret)) && error("ERROR, NAN fount at: ", simparams)

            return ret
        end

        gdmodel = UtilsJL.SimulationUtils.grad_desc_vec(gd_up_fun!;
            smooth = 0.05,
            target, x0, x1, maxΔx, minΔx, gdth, maxsimparams, 
            verbose = true
        )
        z_beta, ug_beta = UtilsJL.SimulationUtils.gd_value(gdmodel)

        # save
        dat = (;PME, z_beta, ug_beta, z_avPME, ug_avPME)
        Dyn.save_maxent(dat, simid, MEmode, simparams)

        break
    end
    
end