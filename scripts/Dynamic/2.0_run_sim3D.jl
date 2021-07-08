@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        SimD3, run_simD3!, check_stst, hist, 
        n_ave_conv!, tail_ave,
        plotsdir, lglob, sglob, Container, vec!, Vi, 
        normalize!, save_sim,
        set_status, get_status, 
        load_batch, save_batch,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS, 
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS
        
    import UtilsJL
    const PltU = UtilsJL.PlotsUtils
    const Ass = UtilsJL.ProjAssistant
    const GU = UtilsJL.GeneralUtils

    using Plots
    using Base.Threads
    using ProgressMeter
    using BenchmarkTools
    using ExtractMacro

    Ass.set_verbose(false)
end

## ------------------------------------------------------
let
    Ds = range(0.1, 0.5; length = 50)
    ϵs = range(0.01, 1.0; length = 50)
    cgs = [15.0, Inf]
    SimD3Id = "SimD3"
    Dyn.sglob(;Ds, ϵs, SimD3Id)
        
    lk = ReentrantLock()
    iter = collect(Iterators.product(Ds, ϵs, cgs))
    @threads for (D, ϵ, cg) in iter

        S = SimD3(;
            # Space
            V = lglob(:Vcell3D), 
            # Chemostat
            D, ϵ, cg,
            X = 0.5, 
            sg = cg,
            Δt = 0.05,
            niters = 15000
        )

        V = S.V 
        P = S.P
        Vz = Vi(V, :z)
        Vug = Vi(V, :ug)
        Vuo = Vi(V, :uo)

        simparams = (;S.D, S.ϵ, S.cg)
        
        # status
        status = get_status(SimD3Id, simparams)
        (status != UNDONE_SIM_STATUS) && continue

        # stst params
        stst_w = 1000
        stst_th = 0.01
        # exploded th
        exploded_Xth = 1e4
        dead_Xth = 1e-4

        # time series
        push_frec = 10
        save_frec = 100 * push_frec
        Xts = Container{Float64}()
        sgts = Container{Float64}()
        z_avts = Container{Float64}()
        ug_avts = Container{Float64}()
        uo_avts = Container{Float64}()
        cgD_Xts = Container{Float64}()

        Pzts = Container{Dict{Float64, Float64}}()
        Pugts = Container{Dict{Float64, Float64}}()
        Puots = Container{Dict{Float64, Float64}}()

        tseries = [Xts, sgts, z_avts, ug_avts, uo_avts, cgD_Xts, Pzts, Pugts, Puots]

        # Info
        info_frec = 50
        _time = time()
        eltime = 0.0

        function dobreak()
            
            # stst
            checktime = (S.it > 2 * stst_w) && iszero(rem(S.it, stst_w))
            checktime |= (S.it == S.niters)
            if checktime && check_stst(stst_w, stst_th, vec(Xts))
                status = set_status(STST_SIM_STATUS, SimD3Id, simparams)
                return true
            end

            # dead
            if (S.X < dead_Xth) || 
                    ((S.it == S.niters) && isinf(S.cg) && tail_ave(S.z_avs, 50) < S.D)
                status = set_status(DEAD_SIM_STATUS, SimD3Id, simparams)
                return true
            end
            # explosion
            if (S.X > exploded_Xth) || 
                    ((S.it == S.niters) && isinf(S.cg) && tail_ave(S.z_avs, 50) > S.D)
                status = set_status(EXPLODED_SIM_STATUS, SimD3Id, simparams)
                return true
            end
            
            # ended
            if (S.it == S.niters)
                status = set_status(NITERS_SIM_STATUS, SimD3Id, simparams)
                return true
            end

        end

        function print_info(msg::String)
            @extract S: it cg D ϵ z_av X dXdt sg ug_av cgD_X

            # Info
            thid = threadid()
            zdiff = abs(z_av - D)
            ugdiff = abs(ug_av - cgD_X)
            lock(lk) do
                @info(msg, 
                    it, thid,
                    status,
                    cg, ϵ,
                    (D, z_av, zdiff), 
                    (ug_av, cgD_X, ugdiff),
                    X, dXdt, 
                    sg, 
                    tfactor, last_err,
                    eltime
                )
            end
        end

        function feedback()
            @extract S: it cg D ϵ z_av ug_av uo_av X dXdt sg cgD_X
 
            # info
            eltime = (time() - _time)
            _time = time()
            iszero(rem(it, info_frec)) && print_info("At iter")

            # collecting
            if (it == 1) || iszero(rem(it, push_frec)) || (S.it == S.niters)

                push!(Xts, X)
                push!(z_avts, z_av)
                push!(ug_avts, ug_av)
                push!(uo_avts, uo_av)
                push!(cgD_Xts, cgD_X)
                push!(sgts, sg)

                push!(Pzts, hist(Vz, P))
                push!(Pugts, hist(Vug, P))
                push!(Puots, hist(Vuo, P))

                # save batch
                if iszero(rem(it, save_frec)) || (S.it == S.niters)
                    resize!.(tseries)
                    batch = (; 
                        push_frec, 
                        Xts, sgts, z_avts, ug_avts, 
                        uo_avts, cgD_Xts, Pzts, Pugts, Puots
                    )
                    save_batch(batch, SimD3Id, simparams, S.it)
                    empty!.(tseries)
                end

            end
        end

        # ------------------------------------------------------
        # tranformation
        Uug = maximum(Vug) * 1.1
        tfactor = 1.0
        min_tfactor = 1.0001
        tfactor_sign = -1.0
        tfactor_step = 0.5
        last_err = 0.0
        function tranformation() 
            @extract S: it cg D ug_av X cgD_X
            
            P .= P .* (tfactor .- (Vug ./ Uug))
            normalize!(P)
            ug_av = sum(Vug .* P)
            
            # check tolerance
            err = abs(ug_av - cgD_X) / cgD_X
            
            # move if too much or too few
            (err > last_err) && (tfactor_sign *= -1.0)
            tfactor += tfactor_sign * tfactor_step * err
            tfactor = max(min_tfactor, tfactor)

            # check if not working
            iszero(rem(S.it, 100)) && (err > 0.03) ? (tfactor_step *= 1.5) : (tfactor_step *= 0.5)

            last_err = err
        end

        # ------------------------------------------------------
        # run sim
        run_simD3!(S;
            dobreak,
            tranformation,
            feedback
        )
        
        # ------------------------------------------------------
        # save sim
        save_sim(S, simid, simparams)

        # ------------------------------------------------------
        # info
        print_info("At end")
        
    end # for (D, ϵ)

end

