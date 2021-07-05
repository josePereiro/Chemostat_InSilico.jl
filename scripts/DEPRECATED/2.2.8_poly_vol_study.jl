## ----------------------------------------------------------------------------
function closest(v, a)
    ci = firstindex(a)
    for i in eachindex(a)
        if abs(v - a[i]) < abs(v - a[ci])
            ci = i
        end
    end
    return ci
end

## ----------------------------------------------------------------------------
let
    dyn_Ds = []
    dyn_cgD_Xs = []
    color = []
    dyn_pol_vols = []
    dyn_ϵs = []
    
    ## ------------------------------------------------------------
    # collecting
    for Vl in INDEX[:Vls], τ in INDEX[:τs]
        for ϵ in INDEX[:ϵs], D in INDEX[:Ds]

            status = idxdat([:status], Vl, D, ϵ, τ)
            (status != :stst) && continue
            M = idxdat([:M], Vl, D, ϵ, τ)

            cgD_X = M.cg * M.D / M.X
            cgD_X > 1.0 && continue
            M.sg > 0.1 && continue

            @assert M.D == D
            @info("Collecting", status, (Vl, D, ϵ, τ), (D, cgD_X))
            
            push!(dyn_Ds, D); push!(dyn_cgD_Xs, cgD_X)
            push!(dyn_pol_vols, Dyn.poly_vol(M, D, cgD_X))
            push!(dyn_ϵs, ϵ); push!(color, ϵs_colors[ϵ])
        end
    end; println()

    ## ------------------------------------------------------------
    # Dyn plot
    let
        p = scatter(dyn_pol_vols, dyn_Ds; 
            color, label = "", m = 8, 
            xlabel = "poly vol", ylabel = "D"
        )
        mysavefig(p, "eps_vs_poly_vol")

        p = scatter(dyn_Ds, dyn_pol_vols; 
            color, label = "", m = 8, 
            xlabel = "D", ylabel = "poly vol"
        )
        mysavefig(p, "D_vs_poly_vol")
    end

    ## ------------------------------------------------------------
    # vol board
    mfactor = 1.1
    bins = 100
    M = INDEX[:M0]
    Ds = range(0.001, maximum(dyn_Ds) * mfactor; length = bins)
    cgD_Xs = range(0.01, maximum(dyn_cgD_Xs) * mfactor; length = bins)
    
    pvol_board = Dyn.poly_vol_board(M, Ds, cgD_Xs; bins)

    ## ------------------------------------------------------------
    # Plots

    ## ------------------------------------------------------------
    # board
    let
        p = heatmap(Ds, cgD_Xs, pvol_board'; 
            title = "Polytope volume", label = "", 
            xlabel = "D", ylabel = "cgD_X"
        )
        scatter!(p, dyn_Ds, dyn_cgD_Xs; 
            ms = 8, color, label = ""
        )
        
        mysavefig(p, "dyn_stst_plo_volumen"; bins)
    end

    ## ------------------------------------------------------------
    # invert form
    let
        ipsize = similar(pvol_board)
        for (i, _) in enumerate(eachcol(pvol_board))
            ipsize[:, end - (i - 1)] .= pvol_board[:, i]
        end

        yticks = UJL.get_ticks(reverse(cgD_Xs); l = 5) do cgD_X
            -round((maximum(cgD_Xs) - cgD_X); digits = 3)
        end
        
        p = heatmap(Ds, cgD_Xs, ipsize'; 
            title = "Polytope volume", label = "", 
            xlabel = "D", ylabel = "cgD_X", 
            yticks
        )

        scatter!(p, dyn_Ds, maximum(cgD_Xs) .- dyn_cgD_Xs; 
            ms = 8, color , label = ""
        )
        
        mysavefig(p, "inv_dyn_stst_plo_volumen"; bins)
    end

end

