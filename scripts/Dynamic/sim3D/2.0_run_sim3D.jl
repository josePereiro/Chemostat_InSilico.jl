using ProjAssistant
@quickactivate

# ------------------------------------------------------
@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        run_simD3!, SimD3, hist, 
        n_ave_conv!, tail_ave,
        Container, vec!, Vi, 
        normalizeP!, batch_file,
        set_status, get_status,
        save_batch, save_bfiles, save_simdat,
        UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
        EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
        STST_SIM_STATUS, reset!, 
        STST, is_steady

    import SimTools
    import SimTools: grad_desc, gd_value

    using Base.Threads

end

## ------------------------------------------------------
function run_sim3D(simid, Ds, ϵs, cg; 
        dosave_batches = true, 
        dosave_simdat = true, 
        docollect_Ps = true
    )
    
    iter = Iterators.product(Ds, ϵs, cg)
    tot_sims = length(iter)
    Ch = Channel() do Ch_
        for (simi, (D, ϵ, cg)) in enumerate(iter)
            put!(Ch_, (simi, (D, ϵ, cg)))
        end
    end

    # Space
    V = lglob(Dyn, :Vcell3D)
    Vz = Vi(V, :z)
    Vug = Vi(V, :ug)
    Vuo = Vi(V, :uo)

    @threads for _ in 1:nthreads()
        thid = threadid()

        for (simi, (D, ϵ, cg)) in Ch
            
            # ---------------------------------------------------------------
            # sim globals
            # SimD3
            S = SimD3(;
                # Space
                V,
                # Chemostat
                D, ϵ, cg,
                X = 0.5, 
                sg = cg,
                Δt = 0.05,
                niters = 20000
            )

            P = S.P

            simparams = (;S.D, S.ϵ, S.cg)
            
            # batchs
            batch_files = []

            # Welcome & status
            status = get_status(simid, simparams) 
            @info("At", simid, simparams, thid, status)
            
            (status != UNDONE_SIM_STATUS) && continue

            # exploded th
            z_tail_acc = 0.0
            z_tail_count = 0
            z_tail_len = 20
            exploded_Xth = 1e2
            dead_Xth = 5e-2
            z_D_diff_th = 0.05

            # time series
            push_frec = 10
            save_frec = 100 * push_frec
            
            Xts = Float64[]
            sgts = Float64[]
            z_avts = Float64[]
            ug_avts = Float64[]
            uo_avts = Float64[]
            cgD_Xts = Float64[]

            Pzts = Dict{Float64, Float64}[]
            Pugts = Dict{Float64, Float64}[]
            Puots = Dict{Float64, Float64}[]

            tseries = [Xts, sgts, z_avts, ug_avts, uo_avts, cgD_Xts, Pzts, Pugts]

            # stst params
            stst_w = 300
            stst_th = 0.01
            nbuffs = length((:X, :z, :sg, :ug, :uo))
            stst = STST(nbuffs, stst_w)

            # Info
            info_frec = 20
            _time = time()
            eltime = 0.0

            # ---------------------------------------------------------------
            # dobreak
            isbreak = false
            function dobreak()
                
                # niter
                isniter = (S.it >= S.niters)
                
                # stst
                push!(stst, S.X, S.z_av, S.sg, S.ug_av, S.uo_av)
                dostst1 = isniter
                dostst2 = (S.it >= 2 * stst_w) && iszero(rem(S.it, stst_w))
                dostst = dostst1 || dostst2
                isstst = dostst && is_steady(stst, stst_th)

                # z_tail
                iszero(rem(S.it, z_tail_len)) && 
                    (z_tail_count = 0, z_tail_acc = 0.0)
                z_tail_acc += S.z_av
                z_tail_count += 1
                z_tail = z_tail_acc / z_tail_count
                
                # dead
                isdead1 = (S.X < dead_Xth)
                isdead2 = (isstst || isniter) && ((S.D - z_tail) / S.D) > z_D_diff_th
                isdead = isdead1 || isdead2

                # explosion
                isexplosion1 = (S.X > exploded_Xth)
                isexplosion2 = (isstst || isniter) && ((z_tail - S.D) / S.D) > z_D_diff_th
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
            _round(v) = round(v; sigdigits = 4)
            function print_info(msg::String)
                @extract S: it cg D ϵ z_av X dXdt sg ug_av uo_av cgD_X

                
                # Info
                zdiff = abs(z_av - D)
                ugdiff = abs(ug_av - cgD_X)
                tug_relerr = ugdiff / cgD_X
                D, z_av, zdiff = _round.([D, z_av, zdiff])
                ug_av, uo_av, cgD_X, ugdiff = _round.([ug_av, uo_av, cgD_X, ugdiff])
                @info(msg,
                    (simi, tot_sims),
                    it, thid,
                    status,
                    cg, ϵ,
                    (D, z_av, zdiff), 
                    (ug_av, cgD_X, ugdiff),
                    uo_av, X, dXdt, 
                    sg, isbreak,
                    tfactor, tug_relerr,
                    eltime
                )
            end

            # ---------------------------------------------------------------
            # feedback
            function feedback()
                @extract S: it cg D ϵ z_av ug_av uo_av X dXdt sg cgD_X
    
                # info
                eltime = (time() - _time)
                _time = time()
                showinfo = (it == 1)
                showinfo |= iszero(rem(it, info_frec))
                showinfo && print_info("At iter")

                # collecting
                docollect = isbreak
                docollect |= (it == 1)
                docollect |= iszero(rem(it, push_frec)) 
                docollect &= dosave_batches
                if docollect

                    push!(Xts, X)
                    push!(z_avts, z_av)
                    push!(ug_avts, ug_av)
                    push!(uo_avts, uo_av)
                    push!(cgD_Xts, cgD_X)
                    push!(sgts, sg)

                    if docollect_Ps
                        push!(Pzts, hist(Vz, P))
                        push!(Pugts, hist(Vug, P))
                        push!(Puots, hist(Vuo, P))
                    end

                    # save batch
                    dosave = isbreak
                    dosave |= iszero(rem(it, save_frec))
                    if dosave
                        batchdat = (; 
                            push_frec, 
                            Xts, sgts, z_avts, ug_avts, uo_avts,
                            cgD_Xts, Pzts, Pugts, Puots
                        )

                        bfile_params = (simid, simparams, S.it)
                        bfile = save_batch(batchdat, bfile_params...)
                        push!(batch_files, bfile_params)

                        empty!.(tseries)
                    end
                end
            end

            # ---------------------------------------------------------------
            # P tranformation 
            # globals
            Uug = maximum(Vug) * 1.1
            tkernel = (Vug ./ Uug)
            tauxP = similar(P)
            tfactor = 1.0
            min_tfactor = 1.001
            gdth = 1e-3
            gd_err = 0.0
            
            # gd fun
            function gd_up_fun(gdmodel) 
                @extract S: ug_av cgD_X

                tfactor = gd_value(gdmodel)
                tfactor = max(min_tfactor, tfactor)
                tauxP .= S.P .* (tfactor .- tkernel)
                normalizeP!(tauxP)
                ug_av_ = sum(Vug .* tauxP)

                gd_err = ((cgD_X - ug_av_) / cgD_X)

                return ug_av_
            end

            gd_break_cond(gdmodel) = (gd_err > 0.0 && gd_err < gdth)

            function tranformation() 
                @extract S: simit=it cg ug_av X cgD_X
                isinf(cg) && return

                # set up gradient descend
                tauxP .= P
                maxΔx = 5.0
                x0 = tfactor * 0.9
                x1 = tfactor
                target = cgD_X 
                maxiter = 100
                verbose = false

                gdmodel_ = grad_desc(gd_up_fun; 
                    break_cond = gd_break_cond,
                    target, x0, x1, maxΔx, gdth, maxiter, verbose
                )
                gdit = gdmodel_.iter

                # update P
                S.P .= tauxP

                # info
                showinfo = iszero(rem(simit, info_frec))
                showinfo && @info("At transformation", 
                    thid, simit, gdit, tfactor, gd_err
                )
                
            end

            # ------------------------------------------------------
            # run sim
            run_simD3!(S;
                dobreak,
                tranformation,
                feedback
            )
            
            # ------------------------------------------------------
            # Save last batch
            feedback()

            # ------------------------------------------------------
            # save sim
            dosave_simdat && save_simdat(S, simid, simparams)
            dosave_batches && save_bfiles(batch_files, simid, simparams)

            # ------------------------------------------------------
            # info
            print_info("At end")
        
        end # for (D, ϵ, sg)
    end # for thid
end

## ------------------------------------------------------
# finite cg
let
    Ds = range(0.1, 0.5; length = 12)
    ϵs = [
        10.0.^range(log10(1e-3), log10(0.1); length=6);
        range(0.2, 1.0; length=6);
    ]
    cg = 15.0
    simid = "SimD3"

    params = (;Ds, ϵs, cg)
    sglob(Dyn, params, :SimD3, :params, :finite_cg)
    run_sim3D(simid, Ds, ϵs, cg; 
        docollect_Ps = false
    )
end

# ## ------------------------------------------------------
# # infinite cg
# let
#     Ds = range(0.1, 0.5; length = 50)
#     ϵs = range(0.01, 1.0; length = 50)
#     cg = Inf
#     simid = "SimD3"
    # params = (;Ds, ϵs, cg)
    # sglob(Dyn, params, :SimD3, :params, :infinite_cg)
#     run_sim3D(simid, Ds, ϵs, cg; 
#         dosave_batches = false, 
#         dosave_simdat = false
#     )
# end