# ----------------------------------------------------------------------------
function get_neigs(a, i, j; r)
    M, N = size(a)
    is = max(1, i - r):min(M, i + r)
    js = max(1, j - r):min(N, j + r)
    return @views a[is, js]
end

# ----------------------------------------------------------------------------
function kernel1(a, i, j; r = 3)
    val0 = a[i,j]
    (isnan(val0) || isinf(val0)) && return val0

    neigs = get_neigs(a, i, j; r)
    
    M = -Inf
    for val in neigs
        (isnan(val) || isinf(val)) && continue
        M = max(val, M)
    end
    return M
end

# ----------------------------------------------------------------------------
function kernel2(a, i, j; r = 3)
    val0 = a[i,j]
    (isnan(val0) || isinf(val0)) && return val0

    neigs = get_neigs(a, i, j; r)
    
    sum = 0.0 
    N = 0.0
    for val in neigs
        (isnan(val) || isinf(val)) && continue
        sum += val
        N += 1.0
    end
    return sum / N
end

# ----------------------------------------------------------------------------
@time let
    LP_cache = nothing

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        # rand() > 0.2 && continue

        status = idxdat([:status], Vl, D, ϵ, τ)
        (status != :stst) && continue
        M = idxdat([:M], Vl, D, ϵ, τ)

        cgD_X = M.cg * M.D / M.X
        cgD_X > 1.0 && continue
        M.sg > 0.1 && continue

        bins = 50
        vatp_range, _ = Dyn.vatpvg_ranges(M)
        vg_range, _ = Dyn.vgvatp_ranges(M)

        cell_board0 = Matrix{Float64}(undef, length(vatp_range), length(vg_range))
        for (vatpi, vatp) in enumerate(vatp_range)
            lX = get(M.Xb, vatp, Dict())
            for (vgi, vg) in enumerate(vg_range)
                cell_board0[vatpi, vgi] = log(get(lX, vg, NaN))
            end
        end
        cell_board1 = zero(cell_board0)
        
        # kernel1
        for (vatpi, vatp) in enumerate(vatp_range)
            for (vgi, vg) in enumerate(vg_range)
                k = kernel1(cell_board0, vatpi, vgi; r = 10)
                cell_board1[vatpi, vgi] = k
            end
        end

        # kernel2
        for (vatpi, vatp) in enumerate(vatp_range)
            for (vgi, vg) in enumerate(vg_range)
                k = kernel2(cell_board1, vatpi, vgi; r = 5)
                cell_board0[vatpi, vgi] = k
            end
        end

        p = heatmap(vatp_range, vg_range, cell_board0'; 
            xlabel = "vatp", ylabel = "vg"
        )

         # Dynamic marginal
        δ = 0.08
        LP_cache = isnothing(LP_cache) ? Dyn.vgvatp_cache(M) : LP_cache
        f(vatp, vg) = M.Xb[vatp][vg] / M.X
        DyMs = Dyn.get_marginals(f, M; δ, LP_cache, verbose = false)
        vatp_av = Dyn.av(DyMs["vatp"]) 
        vg_av = Dyn.av(DyMs["gt"])  
        
        hline!(p, [cgD_X]; ls = :dot, color = :black, lw = 2, label = "")
        scatter!(p, [vatp_av], [vg_av]; color = :white, label = "")
        
        mysavefig(p, "cell_density_pol"; Vl, D, ϵ, τ)

        println()
    end
end;

