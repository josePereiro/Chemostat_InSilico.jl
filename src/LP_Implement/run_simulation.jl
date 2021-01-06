# Sim function
function run_simulation!(M::SimModel; 
        on_iter = (it, M) -> false,
        cache_marginf = 50,
        verbose = true,
        verbose_frec = 10,
        force = false,
        force_frac = 100
    )

    # cache
    cache = vgvatp_cache(M; marginf = cache_marginf)

    ## ---------------------------------------------------------
    # extract model constants
    Xb =  M.Xb
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    Vg, Kg, Vl, Kl = M.Vg, M.Kg, M.Vl, M.Kl
    σ, τ, Δt = M.σ, M.τ, M.Δt
    
    ## ---------------------------------------------------------
    # update model
    function update_net!(net) 
        net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * M.sg) / (Kg + M.sg))
        net.ub[vl_idx] = max(net.lb[vl_idx], (Vl * M.sl) / (Kl + M.sl))
    end
    update_net!(M.net) 
    net_pool = map((i) -> deepcopy(M.net), 1:nthreads())
    all_nets = [M.net; net_pool]

    verbose && (prog = Progress(M.niters; desc = "Simulating ... ", dt = 0.5))
    i_vatp_range, vatp_range, vg_ranges = nothing, nothing, nothing
    vg_ub0, vl_ub0 = Inf, Inf
    vatpvgN, vatpN = 0, 0
    recompute_ranges_th = (10.0^(-M.δvg))/3
    wage = similar_board(Xb, 0) # world age
    ittime_prom = 0.0

    for it in 1:M.niters

        ittime = @elapsed begin

            ## ---------------------------------------------------------
            on_iter(it, M) && break # feed back

            ## ---------------------------------------------------------
            # Check if the polytope changed
            politope_changed = abs(M.net.ub[vg_idx] - vg_ub0) > recompute_ranges_th ||
                               abs(M.net.ub[vl_idx] - vl_ub0) > recompute_ranges_th

            ## ---------------------------------------------------------
            lazy_iter = !(politope_changed || force || rem(it, force_frac) == 0)

            ## ---------------------------------------------------------
            # ranges (recompute ranges if the polytope change)
            if !lazy_iter
                i_vatp_range, vatp_range, vg_ranges = ivatpvg_ranges(M, net_pool)
                vatpvgN = sum(length.(values(vg_ranges))) # total regions
                vatpN = length(i_vatp_range)
                vg_ub0, vl_ub0 = M.net.ub[vg_idx], M.net.ub[vl_idx]
            end
            
            ## ---------------------------------------------------------
            # complete boards
            if !lazy_iter
                complete_board!(Xb, i_vatp_range, vg_ranges, 0.0)
                complete_board!(wage, i_vatp_range, vg_ranges, 0)
            end

            ## ---------------------------------------------------------
            # vatp dependent
            z = zeros(vatpN)
            Σvg__X = zeros(vatpN)

            @threads for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]

                lXb = Xb[vatp]
                lwage = wage[vatp]

                Σvg__Xi = 0.0
                for vg in vg_range
                    Σvg__Xi += lXb[vg]
                    lwage[vg] = it
                end
                Σvg__X[vatpi] = Σvg__Xi

                vg0 = first(vg_range)
                z[vatpi] = cache[vatp][vg0][obj_idx]
            end

            Σ_vatp_vg__z_X = sum(z[vatpi] * Σvg__X[vatpi] for (vatpi, vatp) in i_vatp_range)
            term2 = M.ϵ * (Σ_vatp_vg__z_X) / vatpvgN # global average

            @threads for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]
                
                vgN = length(vg_range)
                term1 = (1 - M.ϵ) * (z[vatpi]* Σvg__X[vatpi]) / vgN # local average
                
                lXb = Xb[vatp]
                term1_term2 = term1 + term2

                for vg in vg_range
                    Xᵢ₋₁ = lXb[vg]
                    term3 = Xᵢ₋₁ * (M.D + τ * M.sl + σ) # cellular dead/loose

                    # update X
                    ΔX = (term1_term2 - term3) * Δt
                    Xᵢ = Xᵢ₋₁ + ΔX
                    lXb[vg] = max(0.0, Xᵢ)
                end
            end
            # end
            
            ## ---------------------------------------------------------
            # Update Xb (let unfeasible out)
            # vatp
            if !lazy_iter

                ## ---------------------------------------------------------
                # Compute feasible and unfeasible total X
                Xfea = 0.0
                Xunfea = 0.0
                for (vatp, lXb) in Xb
                    lwage = wage[vatp]
                    for (vg, lX) in lXb
                        if lwage[vg] != it # unfeasible
                            Xunfea += lX
                            lXb[vg] = 0.0
                        else
                            Xfea += lX
                        end
                    end
                end

                ## ---------------------------------------------------------
                # redistribute Xunfea (maintaining the distribution shape)
                # ΔX = (X/Xfea)*Xunfea = X * Xunfea/Xfea
                # Compute unfeasible
                Xfac = Xunfea/ Xfea
                for (vatp, lXb) in Xb
                    lwage = wage[vatp]
                    for (vg, lX) in lXb
                        lXb[vg] += lX * Xfac
                    end
                end
            end #if !lazy_iter
            M.X = sum(sum.(values.(values(Xb)))) # total X

            ## ---------------------------------------------------------
            # concs
            Σ_vatp_vg__vg_X_pool = zeros(nthreads())
            Σ_vatp_vg__vl_X_pool = zeros(nthreads())
            @threads for (vatpi, vatp) in i_vatp_range
                tid = threadid()
                lcache = cache[vatp]
                lXb = Xb[vatp]
                for vg in vg_ranges[vatpi]
                    vl = lcache[vg][vl_idx]
                    lX = lXb[vg]
                    Σ_vatp_vg__vg_X_pool[tid] += vg * lX
                    Σ_vatp_vg__vl_X_pool[tid] += vl * lX
                end
            end
            Σ_vatp_vg__vg_X = sum(Σ_vatp_vg__vg_X_pool)
            Σ_vatp_vg__vl_X = sum(Σ_vatp_vg__vl_X_pool)

            # update sg
            Δsg = (-Σ_vatp_vg__vg_X + M.D * (M.cg - M.sg)) * Δt
            M.sg = max(0.0, M.sg + Δsg)

            # update sl
            Δsl = (-Σ_vatp_vg__vl_X + M.D * (M.cl - M.sl)) * Δt
            M.sl = max(0.0, M.sl + Δsl)

            ## ---------------------------------------------------------
            # Update polytope
            foreach(update_net!, all_nets)

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
                (:force, force),
                (:politope_changed, politope_changed),
                (:lazy_iter, lazy_iter),
                (": ------ ", "------------------"),
                (:X, M.X),
                (:sg, M.sg),
                (:Δsg, Δsg),
                (:sl, M.sl),
                (:Δsl, Δsl),
                (:vg_ub, M.net.ub[vg_idx]),
                (:vl_ub, M.net.ub[vl_idx]),
                (:vatpvgN, vatpvgN),
                (:Xb_len, sum(length.(values(Xb)))),
            ]
        )

    end # for it in 1:niters

    verbose && finish!(prog)

    return M
end