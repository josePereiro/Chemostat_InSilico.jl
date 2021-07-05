## ---------------------------------------------------------------
function plot_join!(p, P::Dict; filter = identity, pltparams...)
    # collect axis
    vatps = collect(keys(P))
    vatplb = maximum(vatps)
    vatpub = maximum(vatps)
    Nvatps = length(vatps)
    sort!(vatps)
    vatp_idxmap = Dict(vatp => i for (i, vatp) in enumerate(vatps))
    # @show vatplb, vatpub Nvatps
    
    vgs = []
    for (vatp, Pvg) in P
        push!(vgs, keys(Pvg)...)
        unique!(vgs)
    end
    sort!(vgs)
    vglb = minimum(vgs)
    vgub = maximum(vgs)
    Nvgs = length(vgs)
    vg_idxmap = Dict(vg => i for (i, vg) in enumerate(vgs))
    # @show vglb, vgub Nvgs

    # makes MAT
    Pmat = fill(NaN, Nvatps, Nvgs)
    maxp = -Inf
    for (vatp, Pvg) in P
        vatpi = vatp_idxmap[vatp]
        for (vg, p) in Pvg
            vgi = vg_idxmap[vg]
            Pmat[vatpi, vgi] = p
            maxp = max(maxp, p)
        end
    end

    Pmat ./= maxp
    Pmat .= filter.(Pmat)

    # Plots
    p = heatmap(vatps, vgs, (Pmat)'; 
        title = "Join distribution", label = "", 
        xlabel = "vatp", ylabel = "vg",
        pltparams...
    )
end
plot_join(P::Dict; filter = identity, pltparams...) = plot_join!(plot(), P; filter, pltparams...)

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
    
    # @threads 
    for thid in 1:nths
        for (Vl, D, ϵ, τ) in Ch

            rand() < 0.95 && continue

            ps = Plots.Plot[]

            # LOAD
            status = dyn_status(Vl, D, ϵ, τ)
            status != :stst && continue

            # Dynamic
            method = :dyn
            JDAT = join_dat(method, Vl, D, ϵ, τ)
            PX = JDAT[:P]
            vg_avPX = JDAT[:vg_av]
            vg_vaPX = JDAT[:vg_va]
            vatp_avPX = JDAT[:vatp_av]
            vatp_vaPX = JDAT[:vatp_va]

            lock(WLOCK) do

                @info(string("Doing...", "-"^50), 
                    method, (Vl, D, ϵ, τ), 
                ); println()
                
                p = plot_join(PX; filter,
                    title = string(method),
                    color = :dense,
                )
                scatter!(p, [vatp_avPX], [vg_avPX]; 
                    xerr = [sqrt(vatp_vaPX)],
                    yerr = [sqrt(vg_vaPX)],
                    label = "", m = 5, color = :black, alpha = 0.7
                )

                # mysavefig(p, "join_heatmap"; method, Vl, D, ϵ, τ)
                push!(ps, deepcopy(p))
            end
            
            # MaxEnt
            for MEmode in ME_MODELS

                JDAT = join_dat(MEmode, Vl, D, ϵ, τ)
                PME = JDAT[:P]
                vg_avPME = JDAT[:vg_av]
                vg_vaPME = JDAT[:vg_va]
                vatp_avPME = JDAT[:vatp_av]
                vatp_vaPME = JDAT[:vatp_va]

                lock(WLOCK) do   
                    
                    @info(string("Doing...", "-"^50), 
                        MEmode, (Vl, D, ϵ, τ), 
                    ); println()
                    
                    p = plot_join(PME; filter,
                        title = string(MEmode),
                        color = :dense,
                    )
                    scatter!(p, [vatp_avPME], [vg_avPME]; 
                        xerr = [sqrt(vatp_vaPME)],
                        yerr = [sqrt(vg_vaPME)],
                        label = "", m = 5, color = :black, alpha = 0.7
                    )

                    # mysavefig(p, "join_heatmap"; MEmode, Vl, D, ϵ, τ)
                    push!(ps, deepcopy(p))
                end
            end # for MEmode

            mysavefig(ps, "join_heatmap"; Vl, D, ϵ, τ)

        end # for Ch
    end # for thid

end