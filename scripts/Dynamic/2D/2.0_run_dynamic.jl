@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        Sim2D, run_simD2!, check_stst, hist, n_ave!,
        plotsdir, lglob, sglob, Container, vec!
        
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

    Ds = range(0.1, 0.5; length = 10)
    ϵs = range(0.01, 1.0; length = 10)
    cgs = [15.0, Inf]
    Dyn.sglob(;Ds, ϵs)

    lk = ReentrantLock()

    iter = collect(Iterators.product(Ds, ϵs, cgs))
    @threads for (D, ϵ, cg) in iter

        S = Sim2D(;
            # Space
            Vcell = lglob(:Vcell2D), 
            # Chemostat
            D, ϵ, cg,
            X = 0.5, 
            sg = cg, 
            Δt = 0.05,
            niters = 15000
        )
        @extract S: V P

        simid = (;S.D, S.ϵ, S.cg)
        simdfile = Dyn.procdir("sim2D", simid, ".jls")

        # stst params
        stst_w = 500
        stst_th = 0.01

        # status
        status = :undone

        # time series
        Xts = Container{Float64}()
        sgts = Container{Float64}()
        ug_avs = Container{Float64}()
        z_avs = Container{Float64}()
        cgD_Xs = Container{Float64}()

        Pzts = Container{Dict{Float64, Float64}}()
        Pugts = Container{Dict{Float64, Float64}}()

        # Info
        info_frec = 50
        _time = time()
        eltime = 0.0

        function dobreak()
            
            # dead or explosion
            (S.X > 1e6) && (status = :exploded; return true)
            (S.X < 1e-6) && (status = :dead; return true)

            # stst
            (S.it > 2 * stst_w) && iszero(rem(S.it, stst_w)) &&
                check_stst(stst_w, stst_th, vec(Xts)) && 
                (status = :stst; return true)

            # ended
            (S.it == S.niters) && (status = :niters; return true)
        end

        function print_info(msg::String)
            @extract S: it cg D ϵ z_av X dXdt sg ug_av

            # Info
            thid = threadid()
            cgD_X = cg * D / X
            
            lock(lk) do
                @info(msg, 
                    it, thid,
                    status,
                    cg, ϵ,
                    (D, z_av), 
                    (ug_av, cgD_X),
                    X, dXdt, 
                    sg, 
                    eltime
                )
            end
        end

        function feedback()
            @extract S: it cg D ϵ z_av X dXdt sg ug_av

            # info
            eltime = (time() - _time)
            _time = time()
            iszero(rem(it, info_frec)) && print_info("At iter")

            # collecting
            cgD_X = cg * D / X
            push!(ug_avs, ug_av)
            push!(cgD_Xs, cgD_X)
            push!(z_avs, z_av)
            push!(Xts, X)
            push!(sgts, sg)

            if it == 1 || iszero(rem(it, 10))
                push!(Pzts, hist(V.z, P))
                push!(Pugts, hist(V.ug, P))
            end
        end

        # ------------------------------------------------------
        # tranformation
        tkerner = 1.01 .- (V.ug ./ V.Uug)
        tranformation() = (P .= P .* tkerner)

        # ------------------------------------------------------
        # run sim
        run_simD2!(S;
            dobreak,
            tranformation,
            feedback
        )
        
        # ------------------------------------------------------
        # info
        print_info("At end")
        
        # ------------------------------------------------------
        # plots
        # post process time series
        n = 100
        ug_avs = n_ave!(vec!(ug_avs), n)
        cgD_Xs = n_ave!(vec!(cgD_Xs), n)
        z_avs = n_ave!(vec!(z_avs), n)
        Xts = n_ave!(vec!(Xts), n)
        sgts = n_ave!(vec!(sgts), n)

        lock(lk) do

            # plots
            prms = (;title = string(status), label = "", xrotation = 45)
            pz = bar(hist(V.z, P); xlabel = "z", ylabel = "proj", prms...)
            vline!(pz, [S.z_av]; color = :blue, prms...)
            vline!(pz, [S.D]; color = :red, prms...)
            
            cgD_X = S.cg * S.D / S.X
            pug = bar(hist(V.ug, P); xlabel = "ug", ylabel = "proj", prms...)
            vline!(pug, [S.ug_av]; color = :blue, prms...)
            vline!(pug, [cgD_X]; color = :red, prms...)

            pug_av = plot(ug_avs; color = :blue, ylabel = "ug_av & cgD_Xs", prms...)
            plot!(pug_av, cgD_Xs; color = :red, prms...)
            
            pz_av = plot(z_avs; ylabel = "z_av", prms...)
            pz_av = hline!(pz_av, [S.D]; color = :red, prms...)
            
            pX = plot(Xts; ylabel = "X", prms...)
            psg = plot(sgts; ylabel = "sg", prms...)
            
            PltU.sfig([pz, pug, pug_av, pX, pz_av, psg], 
                [plotsdir()], "sim2D", simid, ".png"
            )

            dat = (;S, status, ug_avs, cgD_Xs, z_avs, Xts, sgts, Pzts)
            Dyn.sdat(dat, simdfile; verbose = false)

        end # lock(lk) do

    end # for (D, ϵ)

end

## ------------------------------------------------------
# # gifs
    # gifs_ = map(zip(Pzts, Pugts)) do hists
    #     Pz_, Pug_ = hists
    #     pz_ = bar(Pz_; xlabel = "z", ylabel = "proj", prms...)
    #     pug_ = bar(Pug_; xlabel = "ug", ylabel = "proj", prms...)
    #     PltU.sfig([pz_, pug_], tempname(), ".png")
    # end
    # @show typeof(gifs_)
    # @show length(gifs_)
    # PltU.save_gif(gifs_, 
    #     plotsdir("P_evolution", (;S.D, S.ϵ, S.cg, S.Δt), ".gif")
    # )