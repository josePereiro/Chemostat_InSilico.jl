## ----------------------------------------------------------------------------
# Steady state_vs_D Dynamic
let
    f = identity
    fields = [:X, :sg, :sl]
    for Vl in INDEX[:Vls], τ in INDEX[:τs]

        # COLLECT
        dat_dict = Dict()
        for ϵ in INDEX[:ϵs], D in INDEX[:Ds]
            
            status = idxdat([:status], Vl, D, ϵ, τ)
            (status != :stst && status != :death) && continue
            M = idxdat([:M], Vl, D, ϵ, τ)

            @info("Collecting", (Vl, D, ϵ, τ))

            for field in fields
                dat = get!(dat_dict, (ϵ, field), Dict())
                get!(dat, :xs, [])
                get!(dat, :ys, [])
                get!(dat, :colors, [])

                y = getfield(M, field)

                push!(dat[:xs], D)
                push!(dat[:ys], y)
                push!(dat[:colors], Gray(ϵ * 0.8))
            end
        end

        # PLOT
        ps = Dict(
            field => plot(;title = "Dynamic Steady State", 
                xlabel = "D", ylabel = string(field)
            )
            for field in fields
        )
        for ((ϵ, field), dat) in dat_dict
            ylabel = field
            p = ps[field]
            color = dat[:colors]
            plot!(p, dat[:xs], f.(dat[:ys] .+ 1e-8); 
                label = "", lw = 4, alpha = 0.7, color
            )
            scatter!(p, dat[:xs], f.(dat[:ys] .+ 1e-8); 
                label = "", m = 8, alpha = 0.7, color
            )
        end
        
        # SAVING
        for (field, p) in ps
            mysavefig(p, "$(field)_vs_D_vs_ϵ"; Vl, τ)
        end
    end
end
