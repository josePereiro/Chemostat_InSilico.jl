# ---------------------------------------------------------------
function plot_PY!(p, P::Dict; filter = identity, pltparams...)

    PY = Dict{Float64, Float64}()
    for (vatp, Pvg) in P
        vglb = minimum(keys(Pvg))
        maxY = vatp/vglb
        for (vg, p) in Pvg
            Y = vatp/vg
            @assert Y > 0
            @assert !isinf(Y)
            dY = maxY - Y
            @assert dY >= 0 dY

            get!(PY, dY, 0.0)
            PY[dY] += p 
        end
    end

    # normalize
    Z = sum(values(PY))
    for (dY, p) in PY
        PY[dY] = p / Z
    end

    # Plots
    dYs = sort!(collect(keys(PY)))
    ps = [filter(PY[dY]) for dY in dYs]
    histogram!(p, dYs;
        weights = ps,
        xlabel = "Y vatp/vg", 
        ylabel = "prob", 
        label = "",
        bins = 20,
        pltparams...
    )

    av_Y = sum(abs.(values(PY)) .* abs.(keys(PY)))
    vline!(p, [av_Y], ; 
        label = "", lw = 8, color = :black, 
        alpha = 0.7, ls = :dash
    )

    return p
end
plot_PY(P::Dict; filter = identity, pltparams...) = plot_PY!(plot(), P; filter, pltparams...)

## ---------------------------------------------------------------
let
    
    ME_MODELS = [
        ME_FULL_POLYTOPE,
    ]

    stdf = 3.0
    filter(x) = x
    
    # ---------------------------------------------------------------
    # COLLECT DATA
    # store collected data from join distributions
    nths = 3 # working threads
        
    LP_cache = nothing
    WLOCK = ReentrantLock()

    Ch = Channel(nthreads()) do Ch_
        for (Vl, D, ϵ, τ) in EXP_PARAMS
            put!(Ch_, (Vl, D, ϵ, τ))
        end
    end
    
    @threads for thid in 1:nths
        for (Vl, D, ϵ, τ) in Ch

            D < 0.01

            ps = Plots.Plot[]

            # LOAD
            status = dyn_status(Vl, D, ϵ, τ)
            status != :stst && continue

            # Dynamic
            method = :dyn
            JDAT = join_dat(method, Vl, D, ϵ, τ)
            PX = JDAT[:P]

            # ZPX = sum(sum.(values.(values(PX))))
            # av_PX = Dyn.ave_over((vatp_, vg_) -> vatp_, PX)
            # @show av_PX

            lock(WLOCK) do
                # return 

                @info(string("Doing...", "-"^50), 
                    method, (Vl, D, ϵ, τ), 
                ); println()
                
                p = plot_PY(PX; filter,
                    title = string(method)
                )

                # mysavefig(p, "PY_vatp_vg"; method, Vl, D, ϵ, τ)
                push!(ps, deepcopy(p))
            end
            
            # MaxEnt
            for MEmode in ME_MODELS

                JDAT = join_dat(MEmode, Vl, D, ϵ, τ)
                PME = JDAT[:P]

                # ZPM = sum(sum.(values.(values(PME))))
                # @show ZPM
                # av_PME = Dyn.ave_over((vatp_, vg_) -> vatp_, PME)
                # @show av_PME
                # @show abs(av_PME - av_PX)/av_PX
                # println()

                lock(WLOCK) do   
                    # return
                    
                    @info(string("Doing...", "-"^50), 
                        MEmode, (Vl, D, ϵ, τ), 
                    ); println()
                    
                    p = plot_PY(PME; filter,
                        title = string(MEmode),
                    )

                    # mysavefig(p, "PY_vatp_vg"; MEmode, Vl, D, ϵ, τ)
                    push!(ps, deepcopy(p))
                end
            end # for MEmode

            lock(WLOCK) do 
                mysavefig(ps, "PY_vatp_vg"; Vl, D, ϵ, τ)
            end

        end # for Ch
    end # for thid

end