## ----------------------------------------------------------------------------
# Fix D moving eps marginals
let
    Vl = MINDEX[:Vls]|> first  
    τ = MINDEX[:τs] |> first
    Ds = MINDEX[:Ds] |> sort
    Ds = Ds[1:end - 2]
    ϵs = MINDEX[:ϵs] |> sort

    for (fix_dim_name, fix_dim_vals, mov_dim_name, mov_dim_vals, expand) in [ 
        ("D", Ds, "ϵ", ϵs, (fix, mov) -> (Vl, fix, mov, τ)), 
        ("ϵ", ϵs, "D", Ds, (fix, mov) -> (Vl, mov, fix, τ))
    ]

        for fix_dim_val in fix_dim_vals
            ps = Plots.Plot[]
            @time for mov_dim_val in mov_dim_vals
                Vl, D, ϵ, τ = expand(fix_dim_val, mov_dim_val)
                status = MINDEX[:STATUS, Vl, D, ϵ, τ]
                status != :stst && continue

                sparams =(;lw = 4)
                gparams = (;xaxis = nothing, yaxis = nothing, grid = false)
                
                for rxn in ["gt", "biom"]
                    title = string(
                        rxn, " marginal (", 
                        fix_dim_name, ": ",  fix_dim_val, ", ",
                        mov_dim_name, ": ", mov_dim_val, ")"
                    )
                    p = plot(;title, gparams...)
                    for  MODsym in ALL_MODELS
                        plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
                    end
                    push!(ps, p)
                end
                
            end
            isempty(ps) && continue

            layout = (div(length(ps), 2), 2)
            figname = string("fix_", fix_dim_name, "_marginals")
            mysavefig(ps, figname; 
                τ, Vl, Symbol(fix_dim_name) => fix_dim_val, layout
            ); println()
        end # for fix_dim_val
    end
end

## ----------------------------------------------------------------------------
# vatp, vg marginals v2
let

    ALL_MODELS = [
        # ME_Z_OPEN_G_OPEN, 
        # ME_Z_OPEN_G_BOUNDED, 
        # ME_Z_EXPECTED_G_OPEN, 
        # ME_Z_EXPECTED_G_BOUNDED, 
        ME_FULL_POLYTOPE,
        # ME_Z_EXPECTED_G_MOVING,
        # ME_Z_FIXXED_G_OPEN, 
        # ME_Z_FIXXED_G_BOUNDED, 
        # ME_Z_EXPECTED_G_EXPECTED,

        # FBA_Z_OPEN_G_OPEN, 
        # FBA_Z_OPEN_G_BOUNDED, 
        # FBA_Z_FIXXED_G_OPEN, 
        # FBA_Z_FIXXED_G_BOUNDED
    ]

    # LEGEND PLOT
    leg_p = plot(;title = "Legend", xaxis = nothing, yaxis = nothing)
    for  (i, MODsym) in ALL_MODELS |> enumerate
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        plot!(leg_p, fill(length(ALL_MODELS) - i, 3); label = string(MODsym), color, ls, lw = 8)        
    end
    plot!(leg_p, fill(-1, 3); 
        label = "DYNAMIC", color = :black, lw = 8
    )

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        status = MINDEX[:STATUS, Vl, D, ϵ, τ]
        status != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf], lw = 4)
        gparams = (;grid = false)
        for rxn in ["gt", "biom", "resp", "lt"]
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in ALL_MODELS
                plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        M = idxdat([:M0], Vl, D, ϵ, τ)
        
        # POLYTOPE
        p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg")
        Dyn.plot_polborder!(p, M)
        Dyn.plot_poldist!(p, M)
        DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        vatp_av = Dyn.av(DyMs["vatp"]) 
        vg_av = Dyn.av(DyMs["gt"]) 

        kwargs = (;alpha = 0.5, label = "", color = :black)
        hline!(p, [vg_av]; lw = 5, ls = :dash, kwargs...)
        hline!(p, [M.cg * M.D / M.X]; lw = 5, ls = :dot, kwargs...)
        vline!(p, [vatp_av]; lw = 5, ls = :dash, kwargs...)
        push!(ps, p)

        # LEGEND
        push!(ps, leg_p)

        # SAVING
        mysavefig(ps, "Marginals_v2"; Vl, D, ϵ, M.sg, M.sl)

        return
    end
end

## ----------------------------------------------------------------------------
# # selected Marginals
# let

#     Vl = MINDEX[:Vls] |> first
#     τ = MINDEX[:τs] |> first
#     D = (MINDEX[:Ds] |> sort)[5]
#     ϵs = MINDEX[:ϵs] |> sort
#     sparams =(;alpha = 0.8, lw = 10, ylim = [0.0, Inf])
#     gparams = (;grid = false)
    
#     MODELS = [ME_BOUNDED, ME_EXPECTED, ME_CUTTED, ME_Z_OPEN_G_OPEN, FBA_BOUNDED, FBA_OPEN]
#     ps = Plots.Plot[]
#     for ϵ in ϵs
#         for rxn in ["gt", "vatp"]
#             MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
#             @info("Doing", (rxn, Vl, D, ϵ, τ))
#             p = plot(;title = "ϵ = $ϵ", 
#                 xlabel = rxn == "gt" ? "vg" : rxn, 
#                 ylabel = "prob", gparams...
#             )
            
#             @time begin
#                 plot_marginals!(p, MODELS, rxn, Vl, D, ϵ, τ; sparams...)
#             end
#             push!(ps, p)
#         end
#     end
#     layout = 4, 2
#     mysavefig(ps, "dyn_vs_model_marginals"; layout, Vl, D, τ)
# end
