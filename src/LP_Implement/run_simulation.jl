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
    damp, Xmin, Xmax = M.damp, M.Xmin, M.Xmax
    
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
    lastΔXb = Dict{Float64, Dict{Float64, Float64}}()
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

            # vatp dependent
            z = zeros(vatpN)
            Σvg__X = zeros(vatpN)

            for (vatpi, vatp) in i_vatp_range
                vg_range = vg_ranges[vatpi]

                lXb = haskey(Xb, vatp) ? Xb[vatp] : (Xb[vatp] = Dict{Float64, Float64}())
                llastΔXb = haskey(lastΔXb, vatp) ? lastΔXb[vatp] : (lastΔXb[vatp] = Dict{Float64, Float64}())
                
                Σvg__Xi = 0.0
                for vg in vg_range
                    Σvg__Xi += haskey(lXb, vg) ? lXb[vg] : (lXb[vg] = Xmin)
                    haskey(llastΔXb, vg) ? llastΔXb[vg] : (llastΔXb[vg] = 0.0)
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
                    Xᵢ = (damp * Xᵢ₋₁) + (1 - damp) * (Xᵢ₋₁ + (ΔX + lastΔX)/2)
                    lXb[vg] = clamp(Xᵢ, Xmin, Xmax)
                    llastΔXb[vg] = ΔX
                end
            end
            # end
            
            ## ---------------------------------------------------------
            # Update Xb (let unfeasible out)
            # vatp
            if politope_changed
                unfea_vatp = setdiff(keys(Xb), vatp_range)
                foreach((vatp) -> delete!(Xb, vatp), unfea_vatp)
                
                # vg
                for (vatpi, vatp) in i_vatp_range
                    lXb = Xb[vatp]
                    unfea_vg = setdiff(keys(lXb), vg_ranges[vatpi])
                    foreach((vg) -> delete!(lXb, vg), unfea_vg)
                end
            end
            M.X = sum(sum.(values.(values(Xb))))

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
            M.sg = max(0.0, damp * M.sg + (1 - damp) * (M.sg + (Δsg + lastΔsg)/2))
            lastΔsg = Δsg
            # update sl
            Δsl = -Σ_vatp_vg__vl_X + M.D * (M.cl - M.sl)
            M.sl = max(0.0, damp * M.sl + (1 - damp) * (M.sl + (Δsl + lastΔsl)/2))
            lastΔsl = Δsl

            ## ---------------------------------------------------------
            # Update polytope
            update_net!(M.net) 
            foreach(update_net!, net_pool)

        end # t = @elapsed begin

        ## ---------------------------------------------------------
        ittime_prom = (ittime_prom * (it - 1) + ittime)/ it

        ## ---------------------------------------------------------
        # Verbose
        verbose && update!(prog, it; showvalues = [
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