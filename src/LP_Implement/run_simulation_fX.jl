# Sim function
function run_simulation_fX!(M::SimModel; 
        on_iter = (it, M) -> false,
        verbose = true,
        verbose_frec = 10,
        LP_cache = nothing,
        cache_marginf::Float64 = 1.5
    )

    # LP_cache
    if(isnothing(LP_cache))
        LP_cache = vgvatp_cache(M; marginf = cache_marginf)
    end

    ## ---------------------------------------------------------
    # extract model constants
    Xb =  M.Xb
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    cg, cl, Vg, Vl = M.cg, M.cl, M.Vg, M.Vl

    σ, τ, Δt = M.σ, M.τ, M.Δt

    verbose && (prog = Progress(M.niters; desc = "Simulating ... ", dt = 0.5))
    ittime_prom = 0.0

    ## ---------------------------------------------------------
    # ranges (recompute ranges if the polytope change)
    vatp_range, vg_ranges = vatpvg_ranges(M)
    i_vatp_range = vatp_range |> enumerate |> collect
    vatpvgN = sum(length.(values(vg_ranges))) # total regions
    vatpN = length(i_vatp_range)

    for it in 1:M.niters

        ittime = @elapsed begin

            ## ---------------------------------------------------------
            on_iter(it, M) && break # feed back

            ## ---------------------------------------------------------
            # Growth average
            av_zX = 0.0
            for (vatp, lX) in M.Xb
                lcache = LP_cache[vatp]
                for (vg, X) in lX
                    z = lcache[vg][obj_idx]
                    av_zX += z * X
                end
            end
            av_zX /= vatpvgN

            ## ---------------------------------------------------------
            # update biomass distribution
            neg_factor = M.D + τ * M.sl + σ
            for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]

                lXb = M.Xb[vatp]
                lcache = LP_cache[vatp]

                for vg in vg_range
                    z = lcache[vg][obj_idx]
                    Xᵢ₋₁ = lXb[vg]
                    neg_term = Xᵢ₋₁ * neg_factor
                    pos_term = (1 - M.ϵ) * Xᵢ₋₁ * z  + M.ϵ * av_zX
                    ΔX = (pos_term - neg_term) * Δt
                    Xᵢ = Xᵢ₋₁ + ΔX
                    lXb[vg] = max(0.0, Xᵢ)
                end
            end

            ## ---------------------------------------------------------
            # update X
            M.X = sum(sum.(values.(values(Xb)))) # total X

            ## ---------------------------------------------------------
            # Enforce uptake balances
            # The total biomass in the culture must be so av_ui <= ci D / X
            av_vg = 0.0
            av_vl = 0.0
            for (vatp, lX) in M.Xb
                lcache = LP_cache[vatp]
                for (vg, X) in lX
                    vl = lcache[vg][vl_idx]
                    av_vg += vg * X
                    av_vl += vl * X
                end
            end

            # if γ < 1.0 biomass need a rectification
            γvg = M.cg == 0.0 || av_vg <= 0.0 ? Inf : (M.cg / av_vg) * (M.D / M.X)
            γvl = M.cl == 0.0 || av_vl <= 0.0 ? Inf : (M.cl / av_vl) * (M.D / M.X)
            γ = min(γvg, γvl) 
            if 0.0 <= γ < 1.0 
                for (vatpi, vatp) in i_vatp_range
                    vg_range = vg_ranges[vatpi]
                    lX = M.Xb[vatp]
                    for vg in vg_range
                        lX[vg] *= γ
                    end
                end

                # recompute 
                M.X *= γ
                av_vg *= γ
                av_vl *= γ
            end
            
            ## ---------------------------------------------------------
            # update vessel concentration
            # update sg
            Δsg = (-av_vg + M.D * (M.cg - M.sg)) * Δt
            M.sg = max(0.0, M.sg + Δsg)

            # update sl
            Δsl = (-av_vl + M.D * (M.cl - M.sl)) * Δt
            M.sl = max(0.0, M.sl + Δsl)

        end # t = @elapsed begin

        ## ---------------------------------------------------------
        # Verbose
        up = verbose && rem(it, verbose_frec) == 0
        up && update!(prog, it; 
            showvalues = [
                (:it, it),
                (:ittime, ittime),
                (:D, M.D),
                (:Δt, Δt),
                (:γ, γ),
                (:γvg, γvg),
                (:γvl, γvl),
                (": ------ ", "------------------"),
                (:X, M.X),
                (:av_zX, av_zX),
                (": ------ ", "------------------"),
                (:sg, M.sg),
                (:Δsg, Δsg),
                (:av_vg, av_vg),
                (:av_vg_ub, M.cg * M.D / M.X),
                (": ------ ", "------------------"),
                (:sl, M.sl),
                (:Δsl, Δsl),
                (:av_vl, av_vl),
                (:av_vl_ub, M.cl * M.D / M.X),
                (": ------ ", "------------------"),
                (:vg_ub, M.net.ub[vg_idx]),
                (:vl_ub, M.net.ub[vl_idx]),
                (:vatpvgN, vatpvgN)
            ]
        )

    end # for it in 1:niters

    verbose && finish!(prog)

    return M
end