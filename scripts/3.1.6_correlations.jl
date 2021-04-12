## ----------------------------------------------------------------------------
# Steady State Model Dynamic correlation
let
    f(x) = log10(abs(x) + 1e-8)

    # PARAMS
    ϵs = MINDEX[:ϵs]
    sim_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    sim_params = collect(sim_params)[1:5:end]
    
    ps = Plots.Plot[]
    for ϵ in ϵs |> sort

        for  MODsym in ALL_MODELS
            
            p = plot(;title = string(MODsym, " ϵ: ", ϵ), 
                xlabel = "dym flxs", ylabel = "model flxs", 
                legend = :topleft
            )
            
            color = MOD_COLORS[MODsym]
            arr = Vector{Float64}(undef, length(sim_params) * length(Dyn.RXNS))
            DYN_flxs, DYN_errs = arr, copy(arr)
            M_flxs, M_errs = copy(arr), copy(arr)
            
            @info("Doing", ϵ, MODsym)
            
            for (i, (Vl, D, τ)) in sim_params |> enumerate
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                
                for rxn in Dyn.RXNS
                    DYN_flx = Dyn.av(DyMs[rxn])
                    DYN_err = Dyn.va(DyMs[rxn]) |> sqrt
                    ME_flx = Dyn.av(Ms[rxn])
                    ME_err = Dyn.va(Ms[rxn]) |> sqrt
                    (isnan(DYN_flx) || isnan(ME_flx)) && continue

                    push!(DYN_flxs, DYN_flx)
                    push!(DYN_errs, DYN_err)
                    push!(M_flxs, ME_flx)
                    push!(M_errs, ME_err)
                end
            end
            
            xs = DYN_flxs
            ys = M_flxs
            l = minimum(f.([xs; ys]))            
            u = maximum(f.([xs; ys]))    
            m = abs(u - l) * 0.1        
            
            scatter!(p, f.(xs), f.(ys); ms = 8, alpha = 0.5, color, label = "")
            plot!(p, [l - m, u + m], [l - m, u + m]; label = "", ls = :dash, alpha = 0.8)
            push!(ps, p)
        end # for  MODsym
    end # for ϵ

    # saving
    rows = Int(round(length(ps) / length(ALL_MODELS), RoundUp))
    cols = length(ALL_MODELS)
    layout = rows, cols
    mysavefig(ps, "flxs_corr"; layout) 
end

## ----------------------------------------------------------------------------
# flx corrs
let

    # SETUP
    FLXS = ["vatp", "gt"]
    ϵs = MINDEX[:ϵs]
    color_pool = Dict(
        ϵ => c for (ϵ, c) in 
        zip(ϵs, Plots.distinguishable_colors(length(ϵs)))
    )

    # COLLECT
    dat_pool = Dict()
    for (Vl, D, ϵ, τ) in EXP_PARAMS

        # LOAD
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        @info("Collecting", (Vl, D, ϵ, τ))
        
        for flx in FLXS
            
            dym_vatp = Dyn.av(DyMs[flx])
            # dym_vatp = rand() # Test

            for MODsym in ALL_MODELS
                dat = get!(dat_pool, (flx, MODsym), Dict())
                get!(dat, :xs, [])
                get!(dat, :ys, [])
                get!(dat, :colors, [])

                Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                m_vatp = Dyn.av(Ms[flx])
                # m_vatp = rand() # Test

                push!(dat[:xs], dym_vatp)
                push!(dat[:ys], m_vatp)
                push!(dat[:colors], color_pool[ϵ])
            end
        end
    end

    # PLOT CORRS
    ps_pool = Dict()
    sdat_pool = sort(collect(dat_pool), by = (p) -> last(first(p))) # sort by MODsym
    for ((flx, MODsym), dat) in sdat_pool
        ps = get!(ps_pool, flx, Plots.Plot[])
        p = plot(; title = string("dynamic stst: ", MODsym), 
            xlabel = "dyn $flx", ylabel = "model $flx"
        ) 

        color = dat[:colors]
        scatter!(p, dat[:xs], dat[:ys]; color, 
            label = "", m = 8
        )
        vals = sort([dat[:xs]; dat[:ys]])
        plot!(p, vals, vals; label = "", color = :black, 
            alpha = 0.8, lw = 3, ls = :dash
        )
        push!(ps,p)
        
    end

    # LEGEND
    leg_p = plot(;title = "legend", xaxis = nothing, ylabel = "ϵ")
    for ϵ in ϵs
        color = color_pool[ϵ]
        plot!(leg_p, fill(ϵ, 10); color, label = string("ϵ: ", ϵ), lw = 8)
    end

    # SAVE
    for (flx, ps) in ps_pool
        push!(ps, leg_p)
        mysavefig(ps, "$(flx)_correlation") 
    end
end