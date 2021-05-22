## ---------------------------------------------------------------
function plot_join!(p, P::Dict; filter = identity, pltparams...)
    # collect axis
    vatps = collect(keys(P))
    vatplb = maximum(vatps)
    vatpub = maximum(vatps)
    Nvatps = length(vatps)
    sort!(vatps)
    vatp_idxmap = Dict(vatp => i for (i, vatp) in enumerate(vatps))
    @show vatplb, vatpub Nvatps
    
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
    @show vglb, vgub Nvgs

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
    filter(x) = log2(x)
    
    # ---------------------------------------------------------------
    # COLLECT DATA
    # store collected data from join distributions
    nths = 3 # working threads
        
    LP_cache = nothing
    WLOCK = ReentrantLock()

    Ch = Channel(nthreads()) do Ch_
        for (Vl, D, ϵ, τ) in rand(collect(EXP_PARAMS), 3)
            put!(Ch_, (Vl, D, ϵ, τ))
        end
    end
    
    # @threads 
    for thid in 1:nths
        for (Vl, D, ϵ, τ) in Ch

            rand() < 0.95 && continue

            ps = Plots.Plot[]

            # LOAD
            status = MINDEX[:STATUS, Vl, D, ϵ, τ]
            status != :stst && continue

            # Dynamic
            M0 = idxdat([:M0], Vl, D, ϵ, τ)
            lock(WLOCK) do; isnothing(LP_cache) && (LP_cache = Dyn.vgvatp_cache(M0)) end
            z(vatp, vg) = LP_cache[vatp][vg][M0.obj_idx]
            vatp(vatp, vg) = vatp
            vg(vatp, vg) = vg

            fX(vatp, vg) = M0.Xb[vatp][vg] / M0.X
            PX = Dyn.get_join(fX, M0)
            
            vg_avPX = Dyn.ave_over(vg, PX)
            vg_vaPX = Dyn.ave_over((vatp_, vg_) -> (vg_avPX - vg(vatp_, vg_))^2, PX) 
            vatp_avPX = Dyn.ave_over(vatp, PX)
            vatp_vaPX = Dyn.ave_over((vatp_, vg_) -> (vatp_avPX - vatp(vatp_, vg_))^2, PX) 

            lock(WLOCK) do

                method = :dyn
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
                M = idxdat([MEmode, :M], Vl, D, ϵ, τ)

                beta_biom = idxdat([MEmode, :beta_biom], Vl, D, ϵ, τ)
                beta_vg = idxdat([MEmode, :beta_vg], Vl, D, ϵ, τ)
                fME(vatp_, vg_) = exp((beta_biom * z(vatp_, vg_)) + (beta_vg * vg(vatp_, vg_)))
                PME = Dyn.get_join(fME, M)
                
                vg_avPME = Dyn.ave_over(vg, PME)
                vg_vaPME = Dyn.ave_over((vatp_, vg_) -> (vg_avPME - vg_)^2, PME) 
                vatp_avPME = Dyn.ave_over(vatp, PME)
                vatp_vaPME = Dyn.ave_over((vatp_, vg_) -> (vatp_avPME - vatp(vatp_, vg_))^2, PME) 

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