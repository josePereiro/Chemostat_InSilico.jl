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
    net = M.net
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    Vg, Kg, Vl, Kl = M.Vg, M.Kg, M.Vl, M.Kl
    damp, Xmin = M.damp, M.Xmin
    
    ## ---------------------------------------------------------
    # prepare model
    net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * sg) / (Kg + sg))
    net.ub[vl_idx] = max(net.lb[vl_idx], (Vl * sl) / (Kl + sl))

    verbose && (prog = Progress(M.niters; desc = "Simulating ... "))
    for it in 1:M.niters

        ## ---------------------------------------------------------
        at_iter(it, M) # feed back

        ## ---------------------------------------------------------
        # ranges
        vatp_range, vg_ranges = vatpvg_ranges(M::SimModel)
        vatpvgN = sum(length.(values(vg_ranges))) # total regions

        # vatp dependent
        z = Dict{Float64, Float64}()
        Σvg__X = Dict{Float64, Float64}()

        for (i, vatp) in enumerate(vatp_range)
            @inbounds vg_range = vg_ranges[i]

            lXb = get!(Xb, vatp, Dict{Float64, Float64}())
            
            Σvg__X[vatp] = sum(get!(lXb, vg, Xmin) for vg in vg_range)

            vg0 = first(vg_range)
            z[vatp] = cache[vatp][vg0][obj_idx]
            # vgN[vatp] = length(vg_range)
        end

        Σ_vatp_vg__z_X = sum(z[vatp] * Σvg__X[vatp] for vatp in vatp_range)
        term2 = M.ϵ * (Σ_vatp_vg__z_X) / vatpvgN

        for (i, vatp) in enumerate(vatp_range)
            @inbounds vg_range = vg_ranges[i]

            vgN = length(vg_range)
            term1 = (1 - M.ϵ) * (z[vatp]* Σvg__X[vatp]) / vgN
            
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
        
        ## ---------------------------------------------------------
        # Update Xb (let unfeasible out)
        temp_Xb = deepcopy(Xb)
        empty!(Xb)
        for (i, vatp) in enumerate(vatp_range)
            lXb = Xb[vatp] = Dict{Float64, Float64}() 
            for vg in @inbounds vg_ranges[i]
                lXb[vg] = temp_Xb[vatp][vg]
            end
        end
        X = sum(sum.(values.(values(Xb))))

        ## ---------------------------------------------------------
        # concs
        Σ_vatp_vg__vg_X = 0.0
        Σ_vatp_vg__vl_X = 0.0
        for (i, vatp) in enumerate(vatp_range)
            lcache = cache[vatp]
            lXb = Xb[vatp]
            for vg in @inbounds vg_ranges[i]
                vl = lcache[vg][vl_idx]
                lX = lXb[vg]
                Σ_vatp_vg__vg_X += vg * lX
                Σ_vatp_vg__vl_X += vl * lX
            end
        end
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
        net.ub[vg_idx] = max(net.lb[vg_idx], (M.Vg * sg) / (M.Kg + sg))

        ## ---------------------------------------------------------
        # Verbose
        verbose && update!(prog, it; showvalues = [
                (:it, it),
                (:D, M.D),
                (:damp, damp),
                (": ------ ", "------------------"),
                (:X, X),
                (:sg, sg),
                (:Δsg, Δsg),
                (:sl, sl),
                (:Δsl, Δsl),
                (:vg_ub, net.ub[vg_idx]),
                (:vatpvgN, vatpvgN),
                (:Xb_len, sum(length.(values(Xb)))),
            ]
        )

    end # for it in 1:niters
    verbose && finish!(prog)

    return M
end