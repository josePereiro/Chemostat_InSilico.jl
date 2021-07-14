@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        SimD2, run_simD2!, check_stst, hist, 
        n_ave_conv!, tail_ave,
        plotsdir, lglob, sglob, Container, vec!, Vi, 
        normalize!, save_sim,
        set_status, get_status,
        load_batch, save_batch,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS, reset!
        
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
    # Ds = range(0.1, 0.5; length = 8)
    # ϵs = range(0.01, 1.0; length = 8)
    # cgs = [15.0, Inf]
    Ds = [2.71e-01]
    cgs = [1.50e+01] 
    ϵs = [1.00e+00]
    simid = "SimD2"
    Dyn.sglob(;Ds, ϵs, cgs, SimD2Id = simid)
    
    iter = Iterators.product(Ds, ϵs, cgs)
    tot_sims = length(iter)
    Ch = Channel() do Ch_
        for (simi, (D, ϵ, cg)) in enumerate(iter)
            put!(Ch_, (simi, (D, ϵ, cg)))
        end
    end

    @threads for _ in 1:nthreads()
        thid = threadid()

        for (simi, (D, ϵ, cg)) in Ch

            # SimD2
            S = SimD2(;
                # Space
                V = lglob(:Vcell2D),
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

            simparams = (;S.D, S.ϵ, S.cg)
            
            # batchs
            batch_files = String[]
            
            # status
            status = get_status(simid, simparams)
            # (status != UNDONE_SIM_STATUS) && continue # Test

            # stst params
            stst_max_X = -1
            stst_X_acc = 0.0
            stst_last_it = 0
            stst_w = 2000
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
            cgD_Xts = Container{Float64}()

            Pzts = Container{Dict{Float64, Float64}}()
            Pugts = Container{Dict{Float64, Float64}}()

            tseries = [Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts]

            # Info
            info_frec = 50
            _time = time()
            eltime = 0.0

            # ---------------------------------------------------------------
            # dobreak
            isbreak = false
            function dobreak()
                
                # niter
                isniter = (S.it >= S.niters)
                
                # z_tave
                dotave = (!isempty(z_avts) && isniter)
                z_tave = dotave ? tail_ave(vec(z_avts), 50) : 0.0
                
                # stst
                stst_X_acc += S.X
                stst_max_X = max(stst_max_X, S.X)
                dostst1 = isniter && ((S.it - stst_last_it) > stst_w / 2)
                dostst2 = (S.it > 2 * stst_w) && iszero(rem(S.it, stst_w))
                dostst = dostst1 || dostst2
                isstst = dostst && let
                    # check X stst
                    av_X_ = stst_X_acc / (S.it - stst_last_it)
                    rel_Xdiff_ = abs(av_X_ - S.X) / stst_max_X
                    isstst_ = rel_Xdiff_ < stst_th

                    @info("At isstst check", thid, S.it, rel_Xdiff_, stst_th, isstst_)

                    # reset
                    stst_X_acc = 0.0
                    stst_last_it = S.it
                    
                    # reset
                    isstst_
                end
                
                # dead
                isdead1 = (S.X < dead_Xth)
                isdead2 = (isniter && isinf(S.cg) && z_tave < S.D)
                isdead = isdead1 || isdead2

                # explosion
                isexplosion1 = (S.X > exploded_Xth)
                isexplosion2 = (isniter && isinf(S.cg) && z_tave > S.D)
                isexplosion = isexplosion1 || isexplosion2
                
                # ended
                if isdead
                    status = set_status(DEAD_SIM_STATUS, simid, simparams)
                elseif isexplosion
                    status = set_status(EXPLODED_SIM_STATUS, simid, simparams)
                elseif isstst
                    status = set_status(STST_SIM_STATUS, simid, simparams)
                elseif isniter
                    status = set_status(NITERS_SIM_STATUS, simid, simparams)
                end
                
                isbreak = isdead || isexplosion || isstst || isniter

                isbreak && @info("At isbreak", thid, status, isdead, isexplosion, isstst, isniter)

                return isbreak

            end

            # ---------------------------------------------------------------
            # print_info
            function print_info(msg::String)
                @extract S: it cg D ϵ z_av X dXdt sg ug_av cgD_X

                # Info
                zdiff = abs(z_av - D)
                ugdiff = abs(ug_av - cgD_X)
                @info(msg,
                    (simi, tot_sims),
                    it, thid,
                    status,
                    cg, ϵ,
                    (D, z_av, zdiff), 
                    (ug_av, cgD_X, ugdiff),
                    X, dXdt, 
                    sg, isbreak,
                    tfactor, last_err,
                    eltime
                )
            end

            # ---------------------------------------------------------------
            # feedback
            function feedback()
                @extract S: it cg D ϵ z_av ug_av X dXdt sg cgD_X
    
                # info
                eltime = (time() - _time)
                _time = time()
                showinfo = iszero(rem(it, info_frec))
                showinfo && print_info("At iter")

                # collecting
                docollect = isbreak
                docollect |= (it == 1)
                docollect |= iszero(rem(it, push_frec)) 
                if docollect
                    push!(Xts, X)
                    push!(z_avts, z_av)
                    push!(ug_avts, ug_av)
                    push!(cgD_Xts, cgD_X)
                    push!(sgts, sg)

                    push!(Pzts, hist(Vz, P))
                    push!(Pugts, hist(Vug, P))

                    # save batch
                    dosave = isbreak
                    dosave |= iszero(rem(it, save_frec))
                    if dosave
                        resize!.(tseries)
                        batch = (; 
                            push_frec, 
                            Xts, sgts, z_avts, ug_avts, 
                            cgD_Xts, Pzts, Pugts
                        )
                        save_batch(batch, simid, simparams, S.it)
                        reset!.(tseries)
                    end

                end
            end

            # ------------------------------------------------------
            # tranformation
            Uug = maximum(Vug) * 1.1
            tkernel = (Vug ./ Uug)
            tfactor = 1.0
            min_tfactor = 1.0001
            tfactor_sign = 1.0
            tfactor_step = 0.5
            last_err = 0.0
            function tranformation() 
                @extract S: it cg D ug_av X cgD_X
                isinf(cg) && return
                
                P .= P .* (tfactor .- tkernel)
                normalize!(P)
                ug_av = sum(Vug .* P)
                
                # check tolerance
                err = abs(ug_av - cgD_X) / cgD_X
                
                # move if too much or too few
                (err > last_err) && (tfactor_sign *= -1.0)
                tfactor += tfactor_sign * tfactor_step * err
                tfactor = max(min_tfactor, tfactor)

                # check if not working
                iszero(rem(S.it, 100)) && 
                    (err > 0.03) ? (tfactor_step *= 1.5) : (tfactor_step *= 0.5)

                last_err = err
            end

            # ------------------------------------------------------
            # run sim
            run_simD2!(S;
                dobreak,
                tranformation,
                feedback
            )
            
            # ------------------------------------------------------
            # Save last batch
            feedback()

            # ------------------------------------------------------
            # save sim
            save_sim(S, simid, simparams)

            # ------------------------------------------------------
            # info
            print_info("At end")
        
        end # for (D, ϵ, sg)
    end # for thid
end