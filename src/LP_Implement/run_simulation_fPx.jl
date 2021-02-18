# Sim function
function run_simulation_fPx!(M::SimModel; 
        on_iter = (it, M) -> false,
        verbose = true,
        verbose_frec = 10,
        it0 = 1,
        LP_cache = nothing,
        cache_marginf::Float64 = 0.1
    )

    # LP_cache
    if(isnothing(LP_cache))
        LP_cache = vgvatp_cache(M; marginf = cache_marginf)
    end

    ## ---------------------------------------------------------
    # extract model constants
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    cg, cl, Vg, Vl = M.cg, M.cl, M.Vg, M.Vl
    σ, τ, Δt = M.σ, M.τ, M.Δt
    
    
    verbose && (prog = Progress(M.niters; desc = "Simulating ... ", dt = 0.5))
    ittime_prom = 0.0
    
    ## ---------------------------------------------------------
    # ranges (recompute ranges if the polytope change)
    vatp_range, vg_ranges_i = vatpvg_ranges(M)
    i_vatp_range = vatp_range |> enumerate |> collect
    vg_ranges = Dict(vatp => vg_ranges_i[vatpi] for (vatpi, vatp) in i_vatp_range)

    vatpvgN = sum(length.(values(vg_ranges_i))) # total regions
    vatpN = length(i_vatp_range)
    
    auxXb = deepcopy(M.Xb)
    for it in it0:M.niters

        ittime = @elapsed begin

            ## ---------------------------------------------------------
            # ALGORITHM
            # 1. Update PX
            # 2. Compute total X
            # 3. Check sum_uX ≤ cD/X
            # 4. If 3

            ## ---------------------------------------------------------
            on_iter(it, M) && break # feed back

            ## ---------------------------------------------------------
            # Growth average
            av_zX_PV = 0.0 # Space average
            for (vatp, lX) in M.Xb
                lcache = LP_cache[vatp]
                for (vg, X) in lX
                    z = lcache[vg][obj_idx]
                    av_zX_PV += z * X
                end
            end
            av_zX_PV /= vatpvgN

            ## ---------------------------------------------------------
            # update biomass distribution
            neg_factor = M.D + τ * M.sl + σ
            for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges_i[vatpi]

                lXb = M.Xb[vatp]
                lcache = LP_cache[vatp]

                for vg in vg_range
                    z = lcache[vg][obj_idx]
                    Xᵢ₋₁ = lXb[vg]
                    neg_term = Xᵢ₋₁ * neg_factor
                    pos_term = (1 - M.ϵ) * Xᵢ₋₁ * z  + M.ϵ * av_zX_PV
                    ΔX = (pos_term - neg_term) * Δt
                    Xᵢ = Xᵢ₋₁ + ΔX
                    lXb[vg] = max(0.0, Xᵢ)
                end
            end

            ## ---------------------------------------------------------
            # update X
            M.X = sum(sum.(values.(values(M.Xb)))) # total X

            ## ---------------------------------------------------------
            # Enforce uptake balances
            # The biomass distribution will be moved so it ensure the balance
            sum_vgX = 0.0 # Space total
            sum_vlX = 0.0 # Space total
            for (vatp, lX) in M.Xb
                lcache = LP_cache[vatp]
                for (vg, X) in lX
                    vl = lcache[vg][vl_idx]
                    sum_vgX += vg * X
                    sum_vlX += vl * X
                end
            end

            # if γ == in / up < 1.0 biomass need a rectification
            sum_vgX_ub = M.cg * M.D
            sum_vlX_ub = M.cl * M.D
            γvg = M.cg == 0.0 || sum_vgX <= 0.0 ? Inf : (sum_vgX_ub / sum_vgX)
            γvl = M.cl == 0.0 || sum_vlX <= 0.0 ? Inf : (sum_vlX_ub / sum_vlX)
            γ = min(γvg, γvl) 
            
            # Test
            # sum_vgX_ub = (2 - γ) * sum_vgX
            # sum_vlX_ub = 0.0
            if 0.0 <= γ < 1.0 
                fill_board!(auxXb, 0.0)
                for (vatpi, vatp) in i_vatp_range
                    vg_range = vg_ranges_i[vatpi]
                    lXb = M.Xb[vatp]

                    new_vatp = discretize(γ * vatp, M.δvatp, vatp_range; 
                        def_mode = RoundDown
                    )
                    new_vg_range = vg_ranges[new_vatp]

                    auxlXb = auxXb[new_vatp]
                    for vg in vg_range

                        new_vg = discretize(γ * vg, M.δvg, new_vg_range; 
                            def_mode = RoundDown
                        )
                        auxlXb[new_vg] += lXb[vg]
                    end
                end
                temp = M.Xb
                M.Xb = auxXb
                auxXb = temp
                
                # recompute 
                sum_vgX = 0.0 # Space total
                sum_vlX = 0.0 # Space total
                for (vatp, lX) in M.Xb
                    lcache = LP_cache[vatp]
                    for (vg, X) in lX
                        vl = lcache[vg][vl_idx]
                        sum_vgX += vg * X
                        sum_vlX += vl * X
                    end
                end
            end
            @assert sum_vgX <= sum_vgX_ub
            @assert sum_vlX <= sum_vlX_ub

            ## ---------------------------------------------------------
            # update vessel concentration
            # update sg
            Δsg = (-sum_vgX + M.D * (M.cg - M.sg)) * Δt
            M.sg = max(0.0, M.sg + Δsg)

            # update sl
            Δsl = (-sum_vlX + M.D * (M.cl - M.sl)) * Δt
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
                (:av_zX_PV, av_zX_PV),
                (": ------ ", "------------------"),
                (:sg, M.sg),
                (:Δsg, Δsg),
                (:sum_vgX, sum_vgX),
                (:sum_vgX_ub, sum_vgX_ub),
                (": ------ ", "------------------"),
                (:sl, M.sl),
                (:Δsl, Δsl),
                (:sum_vlX, sum_vlX),
                (:sum_vlX_ub, sum_vlX_ub),
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

## ---------------------------------------------------------
function discretize(v, δ, vrange; def_mode = RoundDown)
    qv = discretize(v, δ; mode = def_mode)
    vL, vU = first(vrange), last(vrange)
    qv < vL && return vL
    qv > vU && return vU
    return qv
end