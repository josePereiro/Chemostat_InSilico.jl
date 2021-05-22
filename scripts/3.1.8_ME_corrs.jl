## ----------------------------------------------------------------------------
# flx corrs
let

    # SETUP
    FLXS = ["vatp", "gt"]
    ϵs = MINDEX[:ϵs]
    Vls = MINDEX[:Vls] |> first
    Ds = MINDEX[:Ds]
    τs = MINDEX[:τs] |> first

    fontsize = 13
    sparams = (;
        dpi = 1000,
        thickness_scaling = 1.3, 
        xguidefontsize = fontsize, 
        yguidefontsize = fontsize
    )
    

    MODELS = [
        ME_FULL_POLYTOPE,
        ME_Z_EXPECTED_G_BOUNDED,
        ME_Z_FIXXED_G_BOUNDED, 
        # ME_Z_OPEN_G_OPEN, 
    ]
    
    # COLLECT
    av_ps_pool = Dict()
    va_ps_pool = Dict()

    for (Vl, D, ϵ, τ) in Iterators.product(Vls, Ds, ϵs, τs)

        # dyn        
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        MOD1 = :dyn
        Ms1 = idxdat([:DyMs], Vl, D, ϵ, τ)
        
        @info("Collecting", (Vl, D, ϵ, τ))

        for MOD2 in MODELS

            Ms2 = idxdat([MOD2, :Ms], Vl, D, ϵ, τ)

            for flx in FLXS

                p0 = plot(;title = "\$ $flx \$", 
                    xlabel = MODEL_LABELS[MOD1], 
                    ylabel = MODEL_LABELS[MOD2]
                ) 
                av_p = get!(av_ps_pool, (MOD1, MOD2, flx), deepcopy(p0))
                va_p = get!(va_ps_pool, (MOD1, MOD2, flx), deepcopy(p0))

                color = ES_COLORS[ϵ]

                # av
                av1 = Dyn.av(Ms1[flx])
                va1 = sqrt(Dyn.va(Ms1[flx]))
                av2 = Dyn.av(Ms2[flx])
                va2 = sqrt(Dyn.va(Ms2[flx]))

                for (p, v1, v2) in [(av_p, av1, av2), (va_p, va1, va2)]
                    scatter!(p, [v1], [v2]; 
                        label = "", color, m = 6, alpha = 0.7,
                        sparams...
                    )
                end

            end
        end
    end

    # display in two columns
    function sortby(pair)
        _, MOD2, flx = first(pair)
        return string(findfirst(isequal(MOD2), MODELS), flx)
    end

    for (fun, pool, limf) in [
        ("av", av_ps_pool, (flx) -> (flx == "gt") ? [0.0, 0.40] : [4.5, 17.5]), 
        ("va", va_ps_pool, (flx) -> (flx == "gt") ? [0.0, 0.13] : [0.0, 5.0])
    ]
        spool = sort!(collect(pool); by = sortby)
        ps = Plots.Plot[]
        for  ((MOD1, MOD2, flx), p) in spool
            plot!(p, x -> x, limf(flx)...; 
                lw = 1, alpha = 0.5, color = :black, 
                label = ""
            )
            push!(ps, p)
        end
        mysavefig(ps, "corrs"; fun)
    end

end