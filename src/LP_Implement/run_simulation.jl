function similar_board(Xb, v0)
    T = typeof(v0)
    sXb = Dict{Float64, Dict{Float64, T}}()
    for (vatp, lX) in Xb
        slX = get!(sXb, vatp, Dict{Float64, T}())
        for (vg, _) in lX
            slX[vg] = v0
        end
    end
    sXb
end

function complete_board!(Xb, i_vatp_range, vg_ranges, deflt)
    T = eltype(values(Xb))
    for (vatpi, vatp) in i_vatp_range
        vg_range = vg_ranges[vatpi]
        lXb = haskey(Xb, vatp) ? Xb[vatp] : (Xb[vatp] = T())
        for vg in vg_range
            !haskey(lXb, vg) && (lXb[vg] = deflt)
        end
    end
    Xb
end



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
    Xb =  M.Xb
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    Vg, Kg, Vl, Kl = M.Vg, M.Kg, M.Vl, M.Kl
    damp, num_min, num_max = M.damp, M.num_min, M.num_max
    
    ## ---------------------------------------------------------
    # update model
    function update_net!(net) 
        net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * M.sg) / (Kg + M.sg))
        net.ub[vl_idx] = max(net.lb[vl_idx], (Vl * M.sl) / (Kl + M.sl))
    end
    update_net!(M.net) 
    net_pool = map((i) -> deepcopy(M.net), 1:nthreads())

    verbose && (prog = Progress(M.niters; desc = "Simulating ... "))
    i_vatp_range, vatp_range, vg_ranges = nothing, nothing, nothing
    vg_ub0, vl_ub0 = Inf, Inf
    vatpvgN, vatpN = 0.0, 0.0
    recompute_ranges_th = (10.0^(-M.θvg))/3
    lastΔXb = similar_board(Xb, 0.0)
    wage = similar_board(Xb, 0) # world age
    lastΔsg = 0.0
    lastΔsl = 0.0
    ittime_prom = 0.0
    for it in 1:M.niters

        ittime = @elapsed begin

            ## ---------------------------------------------------------
            at_iter(it, M) # feed back

            ## ---------------------------------------------------------
            # Check if the polytope changed
            politope_changed = abs(M.net.ub[vg_idx] - vg_ub0) > recompute_ranges_th ||
                               abs(M.net.ub[vl_idx] - vl_ub0) > recompute_ranges_th

            ## ---------------------------------------------------------
            # ranges (recompute ranges if the polytope change)
            if politope_changed
                i_vatp_range, vatp_range, vg_ranges = ivatpvg_ranges(M, net_pool)
                vatpvgN = sum(length.(values(vg_ranges))) # total regions
                vatpN = length(i_vatp_range)
                vg_ub0, vl_ub0 = M.net.ub[vg_idx], M.net.ub[vl_idx]
            end

            ## ---------------------------------------------------------
            # complete boards
            if politope_changed
                complete_board!(Xb, i_vatp_range, vg_ranges, num_min)
                complete_board!(lastΔXb, i_vatp_range, vg_ranges, 0.0)
                complete_board!(wage, i_vatp_range, vg_ranges, 0)
            end

            @assert sum(length.(values(Xb))) == sum(length.(values(lastΔXb)))
            @assert sum(length.(values(Xb))) == sum(length.(values(wage)))
            
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
            term2 = M.ϵ * (Σ_vatp_vg__z_X) / vatpvgN

            @threads for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]
                
                vgN = length(vg_range)
                term1 = (1 - M.ϵ) * (z[vatpi]* Σvg__X[vatpi]) / vgN
                
                lXb = Xb[vatp]
                llastΔXb = lastΔXb[vatp]
                for vg in vg_range
                    Xᵢ₋₁ = lXb[vg]
                    lastΔX = llastΔXb[vg]
                    term3 = Xᵢ₋₁ * M.D
                    # update X
                    ΔX = term1 + term2 - term3
                    # ΔX = clamp(ΔX, num_min, num_max)
                    Xᵢ = (damp * Xᵢ₋₁) + (1.0 - damp) * (Xᵢ₋₁ + (ΔX + lastΔX)/2)
                    lXb[vg] = clamp(Xᵢ, num_min, num_max)
                    llastΔXb[vg] = ΔX
                end
            end
            # end
            
            ## ---------------------------------------------------------
            # Update Xb (let unfeasible out)
            # vatp
            if politope_changed

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
            end
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
            Δsg = -Σ_vatp_vg__vg_X + M.D * (M.cg - M.sg)
            # Δsg = clamp(Δsg, num_min, num_max)
            M.sg = max(0.0, damp * M.sg + (1.0 - damp) * (M.sg + (Δsg + lastΔsg)/2))
            lastΔsg = Δsg
            # update sl
            Δsl = -Σ_vatp_vg__vl_X + M.D * (M.cl - M.sl)
            # Δsl = clamp(Δsl, num_min, num_max)
            M.sl = max(0.0, damp * M.sl + (1.0 - damp) * (M.sl + (Δsl + lastΔsl)/2))
            lastΔsl = Δsl

            ## ---------------------------------------------------------
            # Update polytope
            foreach(update_net!, [M.net; net_pool])

        end # t = @elapsed begin

        ## ---------------------------------------------------------
        ittime_prom = (ittime_prom * (it - 1) + ittime)/ it

        ## ---------------------------------------------------------
        # Verbose
        verbose && update!(prog, it; 
            showvalues = [
                (:it, it),
                (:ittime, ittime),
                (:ittime_prom, ittime_prom),
                (:D, M.D),
                (:damp, damp),
                (:politope_changed, politope_changed),
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