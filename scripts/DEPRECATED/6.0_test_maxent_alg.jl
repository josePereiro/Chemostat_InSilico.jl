using ProjAssistant
@quickactivate

## ------------------------------------------------------
@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
    import Chemostat_InSilico.Dynamic: 
        hist, 
        Vi, normalizeP!, match_moms, maxent,
        algtest_file, save_algtest, load_algtest
        
    import SimTools
    const ST = SimTools

    using Plots
    using ColorSchemes
    import GR
    !isinteractive() && GR.inline("png")

    using Base.Threads
    using ProgressMeter
    using BenchmarkTools
    using ExtractMacro
end

## ------------------------------------------------------
# Utils
# maxent (;PME, z_beta, ug_beta, z_avPME, ug_avPME)
function collect_dat(simid, Ds, cgD_Xs, dat_id::Symbol)
    dat = []
    for D in Ds
        for cgD_X in cgD_Xs
            simparams = (;D, cgD_X)
            me_out = load_algtest(simid, simparams)
            datum = getfield(me_out, dat_id)
            push!(dat, datum)
        end
    end
    return dat
end

## ------------------------------------------------------
# generate maxent data
let
    # Test
    Dyn.@lglob net=net2D V=Vcell2D simid=SimD2Id

    Lz, Uz = Dyn.fva(net, :z)
    marg_z = (Uz - Lz) * 0.1
    Ds = range(Lz + marg_z, Uz - marg_z; length = 5)
    
    Lug, Uug = Dyn.fva(net, :ug)
    marg_ug = (Uug - Lug) * 0.1
    cgD_Xs = range(Lug + marg_ug, Uug - marg_ug; length = 20)
    Dyn.@sglob algtest.Ds algtest.cgD_Xs

    @threads for D in collect(Ds)
        for cgD_X in cgD_Xs

            isfea = Dyn.is_feasible_stst(net, D, cgD_X)
            !isfea && continue

            thid = threadid()
            @info("At", D, cgD_X, thid, )

            simparams = (;D, cgD_X)
            difle = algtest_file(simid, simparams)
            # isfile(difle) && continue

            # maxent (;PME, z_beta, ug_beta, z_avPME, ug_avPME)
            me_out = match_moms(V, D, cgD_X; 
                # Top loop
                top_maxiter = 50,

                # GD globals
                gd_maxiter = 1000,
                gdth = 1e-2,
                verbose = true,

                # stst
                stst_th = 1e-2,
                stst_w = 10
            )

            save_algtest(me_out, simid, simparams)

        end # for cgD_X
    end # @threads for D 

    exit()
end

## ------------------------------------------------------
# generate beta mats
let
    
    Dyn.@lglob net=net2D V=Vcell2D simid=SimD2Id

    _ifnan(d) = (v) -> isnan(v) ? d : v
    _gamma(v) = (minimum(_ifnan(Inf), v), maximum(_ifnan(-Inf), v))
    z_beta0, z_beta1 = _gamma(collect_dat(simid, Ds, cgD_Xs, :z_beta))
    ug_beta0, ug_beta1 = _gamma(collect_dat(simid, Ds, cgD_Xs, :ug_beta))
    
    @show z_beta0, z_beta1
    @show ug_beta0, ug_beta1

    bins = 3000
    z_marg = (z_beta1 - z_beta0) * 0.1
    z_betas = range(z_beta0 - z_marg, z_beta1 + z_marg; length = bins)
    ug_marg = (ug_beta1 - ug_beta0) * 0.1
    ug_betas = range(ug_beta0 - ug_marg, ug_beta1 + ug_marg; length = bins)

    V = Dyn.lglob(:Vcell2D)
    Vz = Vi(V, :z)
    Vug = Vi(V, :ug)
    
    # to avoid exp to explode
    z0 = sum(Vz) / length(Vz)
    ug0 = sum(Vug) / length(Vug)
    Vz0 = Vz .- z0 
    Vug0 = Vug .- ug0

    PME_pool = [ones(length(V)) for th in 1:nthreads()]
    aux_pool = [ones(length(V)) for th in 1:nthreads()]

    dfile = algtest_file(simid, :beta_mats)
    if isfile(dfile)
        H_mat, z_av_mat, ug_av_mat = Dyn.lprocdat(dfile; verbose = false)
    else
        H_mat = fill(-1.0, length(z_betas), length(ug_betas))
        z_av_mat = zeros(length(z_betas), length(ug_betas))
        ug_av_mat = zeros(length(z_betas), length(ug_betas))
    end

    # sim utils
    lk = ReentrantLock()
    z_beta_c = 0
    save_frec = max(1, div(length(z_betas), 100))
    dat = (;z_betas, ug_betas, H_mat, z_av_mat, ug_av_mat)

    @threads for (z_betai, z_beta) in collect(enumerate(z_betas))

        thid = threadid()
        PME = PME_pool[thid]
        z_aux = aux_pool[thid]
        
        # info
        lock(() -> (z_beta_c += 1), lk)
        @info("At", z_beta_c, z_betai, z_beta, thid)

        # check cache
        (H_mat[z_betai, end] != -1.0) && continue
        
        z_aux .= exp.(z_beta .* Vz0)
        @inbounds for (ug_betai, ug_beta) in enumerate(ug_betas)


            PME .= z_aux .* exp.(ug_beta .* Vug0)
            normalizeP!(PME)

            H, z_av, ug_av = 0.0, 0.0, 0.0
            for i in eachindex(PME)
                p, z, ug = PME[i], Vz[i], Vug[i]
                H -= log(p) * p
                z_av += p * z 
                ug_av += p * ug
            end
            
            H_mat[z_betai, ug_betai] = H
            z_av_mat[z_betai, ug_betai] = z_av
            ug_av_mat[z_betai, ug_betai] = ug_av
        end

        # save cache
        lock(lk) do
            dosave = (z_beta_c > 1) && iszero(rem(z_beta_c, save_frec))
            dosave && save_algtest(dat, simid, :beta_mats)
            dosave && @info("Saved cache", z_beta_c, z_betai, z_beta, thid)
        end
    end

    save_algtest(dat, simid, :beta_mats)
end

## ------------------------------------------------------
# valid mats
let
    # maxent dat
    Dyn.@lglob simid=SimD2Id algtest.Ds algtest.cgD_Xs

    # beta mat dat
    z_betas, ug_betas, H_mat, z_av_mat, ug_av_mat = load_algtest(simid, :beta_mats)
    
    iter = collect(Iterators.product(Ds, cgD_Xs))
    @threads for (D, cgD_X) in iter
        
        thid = threadid()
        
        @info("Doing", D, cgD_X, thid)
        
        val_is = Int[]
        rtol = 0.05
        while isempty(val_is) && rtol < 1.5
            @info("Doing", rtol, thid)
            for i in eachindex(H_mat)
                H, z_av, ug_av = H_mat[i], z_av_mat[i], ug_av_mat[i]
                isvalid = !isnan(H) && !isnan(z_av) && !isnan(ug_av)
                isvalid &= isapprox(z_av, D; rtol)
                isvalid &= ug_av <= cgD_X
                isvalid && push!(val_is, i)
            end
            rtol *= 1.1
        end

        save_algtest(val_is, simid, :valid_mat, (;D, cgD_X))
    end
end

## ------------------------------------------------------
# Plots

## ------------------------------------------------------
# valid beta heat maps
let
    # maxent dat
    Dyn.@lglob simid=SimD2Id algtest.Ds algtest.cgD_Xs
    @show simid, Ds, cgD_Xs
    return

    # beta mat dat
    z_betas, ug_betas, H_mat, z_av_mat, ug_av_mat = load_algtest(simid, :beta_mats)

    # # ------------------------------------------------------------------------
    # # beta mats
    # ps = Plots.Plot[]
    # subbins = 600
    # for (matid, mat) in [("H", H_mat), ("z", z_av_mat), ("ug", ug_av_mat)]
    #     ixds0 = 1:max(1, div(size(mat, 1), subbins)):size(mat, 1)
    #     ixds1 = 1:max(1, div(size(mat, 2), subbins)):size(mat, 2)
    #     p = plot(;xlabel = "z_beta", ylabel = "ug_beta", title = matid)
    #     heatmap!(p, z_betas[ixds0], ug_betas[ixds1], mat[ixds0, ixds1]'; label = "")
    #     push!(ps, p)
    # end

    # ffile = Dyn.sfig(ps, 
    #     simid, "algtest", "beta_mats", ".png"
    # )
    # @info("Done", ffile)

    # ------------------------------------------------------------------------
    # collect corrs data
    direct_Hs, direct_zs, direct_ugs = [], [], []
    maxent_Hs, maxent_zs, maxent_ugs = [], [], []
    direct_z_betas, direct_ug_betas = [], []
    maxent_z_betas, maxent_ug_betas = [], []

    # Temp, To delete
    V = Dyn.lglob(:Vcell2D)
    Vz = Vi(V, :z)
    Vug = Vi(V, :ug)
    
    for D in Ds, cgD_X in cgD_Xs
            
        # load valids
        valid_mat = zeros(size(H_mat)...)
        val_is = load_algtest(simid, :valid_mat, (;D, cgD_X))
        valid_mat[val_is] .= H_mat[val_is]
        
        # mat 
        max_H_z_betai, max_H_ug_betai = last(findmax(valid_mat)).I
        direct_z_beta = z_betas[max_H_z_betai]
        direct_ug_beta = ug_betas[max_H_ug_betai]
        
        # maxent
        me_out = load_algtest(simid, (;D, cgD_X))
        @extract me_out: PME maxent_z_beta=z_beta maxent_ug_beta=ug_beta maxent_z=z_avPME maxent_ug=ug_avPME

        # # ------------------------------------------------------------------------
        # # plot valid map
        # p = plot(;xlabel = "z_beta", ylabel = "ug_beta", title = "valid")
        # heatmap!(p, z_betas, ug_betas, valid_mat'; label = "")
        # scatter!(p, [direct_z_beta], [direct_ug_beta]; label = "", m = (9, :square), c = :blue)
        # scatter!(p, [maxent_z_beta], [maxent_ug_beta]; label = "", m = (9, :star), c = :red)
        
        # ffile = Dyn.sfig(p, 
        #     simid, "valid_H_mat", (;D, cgD_X), ".png"
        # )
        # @info("Done", ffile)

        # collect for corrs
        direct_H = H_mat[max_H_z_betai, max_H_ug_betai]
        direct_z = z_av_mat[max_H_z_betai, max_H_ug_betai]
        direct_ug = ug_av_mat[max_H_z_betai, max_H_ug_betai]

        maxent_H = -sum(PME .* log.(PME))
        
        # Temp, to delete
        maxent_ug = sum(PME .* Vug)
        
        # filter
        (abs(direct_z - maxent_z) / D) > 0.1 && continue
        # (abs(direct_ug - maxent_ug) / cgD_X) > 0.1 && continue

        push!(direct_Hs, direct_H)
        push!(direct_zs, direct_z)
        push!(direct_ugs, direct_ug)
        push!(direct_z_betas, direct_z_beta)
        push!(direct_ug_betas, direct_ug_beta)
        
        push!(maxent_Hs, maxent_H)
        push!(maxent_zs, maxent_z)
        push!(maxent_ugs, maxent_ug)
        push!(maxent_z_betas, maxent_z_beta)
        push!(maxent_ug_betas, maxent_ug_beta)

        @info("Doing", D, cgD_X, (direct_ug, maxent_ug))

        println()

    end # for D, cgD_X

    # ------------------------------------------------------------------------
    # plot corrs
    # corrs
    ps = Plots.Plot[]
    for (datid, direct, maxent) in [
                ("H", direct_Hs, maxent_Hs), 
                ("z", direct_zs, maxent_zs),
                ("ug", direct_ugs, maxent_ugs),
                ("z_beta", direct_z_betas, maxent_z_betas),
                ("ug_beta", direct_ug_betas, maxent_ug_betas),
            ]
        p = plot(;xlabel = "direct", ylabel = "maxent", title = datid)
        scatter!(p, direct, maxent; label = "", m = 5)
        all = sort!([direct; maxent])
        plot!(p, all, all; label = "", 
            color = :black, ls = :dash, lw = 3, alpha = 0.5
        )
        push!(ps, p)
    end

    ffile = Dyn.sfig(ps, 
        simid, "algtest", "beta_corrs", ".png"
    )
    @info("Done", ffile)
end

## ------------------------------------------------------
# ug vs H study
let
    Dyn.@lglob net=net2D algtest.Ds algtest.cgD_Xs
    Dyn.@lglob simid=SimD2Id algtest.Ds algtest.cgD_Xs
    Ds = Ds[1:end - 1]

    maxDs = maximum(Ds)
    @show Ds
    @show maxDs

    colors = cgrad(:blues, max(2, length(Ds)))

    ps = Plots.Plot[]
    for (Di,  D) in enumerate(Ds)

        title = round(D; sigdigits = 3)
        thickness_scaling = 1.6
        ug_av_p = plot(;xlabel = "<ug>", ylabel = "H", title, thickness_scaling)
        z_beta_p = plot(;xlabel = "z_beta", ylabel = "H", title, thickness_scaling)
        ug_beta_p = plot(;xlabel = "ug_beta", ylabel = "H", title, thickness_scaling)

        ugs = []
        Hs = []
        z_betas = []
        ug_betas = []
        for cgD_X in cgD_Xs

            isfea = Dyn.is_feasible_stst(net, D, cgD_X)
            !isfea && continue

            simparams = (;D, cgD_X)

            # me_out (;PME, z_beta, ug_beta, z_avPME, ug_avPME)
            dfile = algtest_file(simid, simparams)
            !isfile(dfile) && continue
            me_out = load_algtest(simid, simparams)

            
            PME = me_out.PME
            H = -sum(log.(PME) .* PME)
            z = me_out.z_avPME
            ug = me_out.ug_avPME
            z_beta = me_out.z_beta
            ug_beta = me_out.ug_beta

            @info("At", simparams, H, z, ug, z_beta, ug_beta)

            # filter
            (abs(z - D) / D) > 0.1 && continue
            (abs(ug - cgD_X) / cgD_X) > 0.1 && continue
            (isnan(H) || isnan(ug)) && continue

            push!(Hs, H)
            push!(ugs, ug)
            push!(z_betas, z_beta)
            push!(ug_betas, ug_beta)
            
        end
        for (p, xs) in [(ug_av_p, ugs), (z_beta_p, z_betas), (ug_beta_p, ug_betas)]
            sids = sortperm(xs)
            plot!(p, xs[sids], Hs[sids]; label = "", color = colors[Di], lw = 3)
            scatter!(p, xs[sids], Hs[sids]; label = "", color = colors[Di], lw = 3)
            push!(ps, p)
        end
    end
    # plot(ug_av_p, z_beta_p, ug_beta_p)
    ffile = Dyn.sfig(ps, 
        simid, "algtest", "H_trends", ".png";
        layout = (length(Ds), 3)
    )
    @info("Done", ffile)
   
end