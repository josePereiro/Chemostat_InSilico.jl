## ----------------------------------------------------------------------------
# _vs_time_vs_ϵ_vs_D
let
    marginf = 0.2
    f = identity
    Ds = INDEX[:Ds]
    fields = [:X_ts, :sg_ts, :sl_ts]
    cparams = (;lw = 4, alpha = 0.7)

    for Vl in INDEX[:Vls], τ in INDEX[:τs]
        ps = Plots.Plot[]
        for D in Ds
            TSs = idxdat([:TS], Vl, D, INDEX[:ϵs], τ)
            Dps = Plots.Plot[]
            for field in fields
                ylabel = replace(string(field), "_ts" => "")
                vals = getfield.(TSs, field) 
                ylim = Dyn.lims(marginf, vals...)
                
                p = plot(;xlabel = "time", ylabel, 
                    title = string("D: ", UJL.sci(D)))
                for (ϵ, val) in zip(INDEX[:ϵs], vals) |> collect |> reverse
                    plot!(p, f.(val); 
                        label = string("ϵ: ", ϵ), color = ϵs_colors[ϵ], cparams...
                    )
                end
                push!(ps, p)
                push!(Dps, p)
            end
            layout = length(fields), 1
            mysavefig(Dps, "time_series_vs_ϵ_vs_D"; layout, D, Vl, τ)
        end
        
        layout = length(fields), length(Ds)
        mysavefig(ps, "time_series_vs_ϵ_vs_D"; layout, Vl, τ)
    end
end
