# Sim function
function run_simulation!(M::SimModel; 
        at_iter = (it, M) -> nothing,
        cache_marginf = 1,
        verbose = true
    )

    # cache
    cache = vgvatp_cache(M; marginf = cache_marginf)

    ## ---------------------------------------------------------
    # extract model
    sl_ts, sg_ts, X_ts, Xb = M.sl_ts, M.sg_ts, M.X_ts, M.Xb
    sg, sl = last(sg_ts), last(sl_ts)
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    Vg, Kg, Vl, Kl = M.Vg, M.Kg, M.Vl, M.Kl
    damp, Xmin = M.damp, M.Xmin
    
    ## ---------------------------------------------------------
    # update model
    function update_net!(net) 
        net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * sg) / (Kg + sg))
        net.ub[vl_idx] = max(net.lb[vl_idx], (Vl * sl) / (Kl + sl))
    end
    update_net!(M.net) 
    net_pool = map((i) -> deepcopy(M.net), 1:nthreads())

    verbose && (prog = Progress(M.niters; desc = "Simulating ... "))
    for it in 1:M.niters

        ittime = @elapsed begin

            ## ---------------------------------------------------------
            at_iter(it, M) # feed back

            ## ---------------------------------------------------------
            # ranges
            i_vatp_range, vatp_range, vg_ranges = ivatpvg_ranges(M, net_pool)
            vatpvgN = sum(length.(values(vg_ranges))) # total regions
            vatpN = length(i_vatp_range)

            # vatp dependent
            z = zeros(vatpN)
            Σvg__X = zeros(vatpN)


            for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]

                lXb = haskey(Xb, vatp) ? Xb[vatp] : (Xb[vatp] = Dict{Float64, Float64}())
                
                Σvg__Xi = 0.0
                for vg in vg_range
                    Σvg__Xi += haskey(lXb, vg) ? lXb[vg] : (lXb[vg] = Xmin)
                end
                Σvg__X[vatpi] = Σvg__Xi

                vg0 = first(vg_range)
                z[vatpi] = cache[vatp][vg0][obj_idx]
            end

            Σ_vatp_vg__z_X = sum(z[vatpi] * Σvg__X[vatpi] for (vatpi, vatp) in i_vatp_range)
            term2 = M.ϵ * (Σ_vatp_vg__z_X) / vatpvgN

            @threads for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]
                
                vgN = length(vg_range)
                term1 = (1 - M.ϵ) * (z[vatpi]* Σvg__X[vatpi]) / vgN
                
                lXb = Xb[vatp]
                for vg in vg_range
                    Xᵢ₋₁ = lXb[vg]
                    term3 = Xᵢ₋₁ * M.D
                    # update X
                    ΔX = term1 + term2 - term3
                    Xᵢ = (damp * Xᵢ₋₁) + (1 - damp) * (Xᵢ₋₁ - ΔX)
                    lXb[vg] = max(Xmin, Xᵢ)
                end
            end
            # end
            
            ## ---------------------------------------------------------
            # Update Xb (let unfeasible out)
            # vatp
            unfea_vatp = setdiff(keys(Xb), vatp_range)
            foreach((vatp) -> delete!(Xb, vatp), unfea_vatp)
            
            # vg
            for (vatpi, vatp) in i_vatp_range
                lXb = Xb[vatp]
                unfea_vg = setdiff(keys(lXb), vg_ranges[vatpi])
                foreach((vg) -> delete!(lXb, vg), unfea_vg)
            end
            X = sum(sum.(values.(values(Xb))))

            ## ---------------------------------------------------------
            # concs
            # Σ_vatp_vg__vg_X = 0.0
            # Σ_vatp_vg__vl_X = 0.0
            Σ_vatp_vg__vg_X_pool = zeros(nthreads())
            Σ_vatp_vg__vl_X_pool = zeros(nthreads())
            @threads for (vatpi, vatp) in i_vatp_range
                tid = threadid()
                lcache = cache[vatp]
                lXb = Xb[vatp]
                for vg in vg_ranges[vatpi]
                    vl = lcache[vg][vl_idx]
                    lX = lXb[vg]
                    # Σ_vatp_vg__vg_X += vg * lX
                    Σ_vatp_vg__vg_X_pool[tid] += vg * lX
                    # Σ_vatp_vg__vl_X += vl * lX
                    Σ_vatp_vg__vl_X_pool[tid] += vl * lX

                end
            end
            Σ_vatp_vg__vg_X = sum(Σ_vatp_vg__vg_X_pool)
            Σ_vatp_vg__vl_X = sum(Σ_vatp_vg__vl_X_pool)

            # update sg
            Δsg = -Σ_vatp_vg__vg_X + M.D * (M.cg - sg)
            sg = max(0.0, damp * sg + (1 - damp) * (sg + Δsg))
            # update sl
            Δsl = -Σ_vatp_vg__vl_X + M.D * (M.cl - sl)
            sl = max(0.0, damp * sl + (1 - damp) * (sl + Δsl))

            ## ---------------------------------------------------------
            # Store
            push!(X_ts, X)
            push!(sl_ts, sl)
            push!(sg_ts, sg)

            ## ---------------------------------------------------------
            # Update polytope
            update_net!(M.net) 
            foreach(update_net!, net_pool)

        end # t = @elapsed begin

        ## ---------------------------------------------------------
        # Verbose
        verbose && update!(prog, it; showvalues = [
                (:it, it),
                (:ittime, ittime),
                (:D, M.D),
                (:damp, damp),
                (": ------ ", "------------------"),
                (:X, X),
                (:sg, sg),
                (:Δsg, Δsg),
                (:sl, sl),
                (:Δsl, Δsl),
                (:vg_ub, M.net.ub[vg_idx]),
                (:vatpvgN, vatpvgN),
                (:Xb_len, sum(length.(values(Xb)))),
            ]
        )

    end # for it in 1:niters

    verbose && finish!(prog)

    return M
end