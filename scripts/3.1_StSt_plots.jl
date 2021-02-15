import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement
        
    using ProgressMeter
    using Plots.Measures

    using Plots
    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    using Base.Threads
    using Dates
    using Statistics
    using InteractiveUtils
    using Random
    using Colors
    using FileIO

end
## ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
minmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

## ----------------------------------------------------------------------------
# Prepare network
const ME_HOMO = :ME_HOMO           # Do not use extra constraints
const ME_EXPECTED = :ME_EXPECTED   # Match ME and Dy biom average
const ME_CUTTED = :ME_CUTTED     # Match ME and Dy biom average and constraint av_ug
const ME_BOUNDED = :ME_BOUNDED     # Fix biom around observed

const FBA_OPEN = :FBA_OPEN
const FBA_BOUNDED = :FBA_BOUNDED

MOD_COLORS = Dict(
    ME_HOMO => :red, 
    ME_BOUNDED => :orange,
    ME_CUTTED => :blue,
    ME_EXPECTED => :purple,
    FBA_BOUNDED => :green,
    FBA_OPEN => :yellow,
)

MOD_LS = Dict(
    ME_HOMO => :dash, 
    ME_BOUNDED => :dash,
    ME_CUTTED => :dash,
    ME_EXPECTED => :dash,
    FBA_BOUNDED => :dot,
    FBA_OPEN => :dot,
)

## ----------------------------------------------------------------------------
# Marginals
# MDAT[MODsym, :M, Vl, D, ϵ, τ]
# MDAT[MODsym, :Ms, Vl, D, ϵ, τ]
# MDAT[MODsym, :beta0, Vl, D, ϵ, τ]
# MDAT[:STATUS, Vl, D, ϵ]
MINDEX = UJL.load_data(InCh.MARGINALS_INDEX_FILE)
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...)
idxdat(dk, indexks...) = InLP.idxdat(MINDEX, dk, indexks...)

## ----------------------------------------------------------------------------
# PLOTS
mysavefig(p, pname; params...) = 
    InLP.mysavefig(p, pname, InLP.DYN_FIGURES_DIR, fileid; params...)

## ----------------------------------------------------------------------------
# Plot functions
function plot_pol!(p, MODsym, Vl, D, ϵ, τ; sparams...)
    
    vatp_range, vg_ranges = idxdat([MODsym, :POL], Vl, D, ϵ, τ)
    vatps, vgLs, vgUs = [], [], []
    
    for (vatpi, vatp) in enumerate(vatp_range)
        vg_range = vg_ranges[vatpi]
        isempty(vg_range) && continue
        vgL, vgU = minmax(vg_range)
        push!(vatps, vatp)
        push!(vgLs, vgL)
        push!(vgUs, vgU)
    end

    color = MOD_COLORS[MODsym]
    params = (;label = "", alpha = 0.8, color, lw = 8, sparams...)
    plot!(p, [vatps], [vgLs]; params...)
    plot!(p, [vatps], [vgUs]; params...)
    p
end

function plot_marginals!(p, MODsyms, rxn, Vl, D, ϵ, τ; sparams...)

    DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
    plot!(p, DyMs[rxn]; label = "", sparams..., color = :black)
    
    # Marginals
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        plot!(p, Ms[rxn]; label = "", color, ls, sparams...)
    end
    return p
end
plot_marginals!(p, MODsyms::Symbol, rxn, Vl, D, ϵ, τ; sparams...) = 
    plot_marginals!(p, [MODsyms], rxn, Vl, D, ϵ, τ; sparams...)


# ## ----------------------------------------------------------------------------
# let
#     Vl = MINDEX[:Vls] |> first
#     τ = MINDEX[:τs] |> first
#     D = MINDEX[:Ds][2]
#     ϵ = MINDEX[:ϵs][end]
#     MODsym = FBA_OPEN

#     @info("Doing", 
#         (Vl, D, ϵ, τ)
#     ); println()
#     M = idxdat([MODsym, :M], Vl, D, ϵ, τ)

#     p = plot()
#     sparams = (;)
#     # plot_polgrid!(p, M; sparams...)
#     InLP.plot_polborder!(p, M; sparams...)
#     InLP.plot_poldist!(p, M; rand_th = 1.0, 
#         hits_count = Int(1e3), maxiters = Int(1e5), 
#         static_th = 0.05,
#         skwargs = (;color = :red, alpha = 0.8)
#     )

#     mysavefig(p, "test")
# end

## ----------------------------------------------------------------------------
# Error
let
    Vl = MINDEX[:Vls] |> first
    τ = MINDEX[:τs] |> first
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort

    MODELS = [ME_BOUNDED, ME_CUTTED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]

    p = plot(;tile = "Prediction error", xlabel = "log ϵ", ylabel = "maximum err")
    for MODsym in MODELS
        xs, ys, yerrs = [], [], []
        for ϵ in ϵs

            errs = []
            @info("Doing", ϵ, MODsym); println()
            for D in Ds
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                
                for rxn in InLP.RXNS
                    DYN_flx = InLP.av(DyMs[rxn])
                    ME_flx = InLP.av(Ms[rxn])
                    (isnan(DYN_flx) || isnan(ME_flx)) && continue

                    err = ((DYN_flx - ME_flx)^2)/abs(DYN_flx)
                    push!(errs, err)
                end
            end

            push!(xs, ϵ)
            push!(ys, maximum(errs))
            # push!(yerrs, std(errs))

        end # for ϵ in ϵs

        noise = xs .* 0.1 .* rand.()
        params = (;alpha = 0.8, color = MOD_COLORS[MODsym])
        plot!(p, log10.(xs .+ noise), ys; 
            label = "", lw = 3, ls = :dash, params...
        )
        scatter!(p, log10.(xs .+ noise), ys; # yerr = yerrs, 
            ms = 8, params..., label = string(MODsym), legend = :topleft
        )
    end #  for MODsym 

    mysavefig(p, "eps_vs_err")

end

## ----------------------------------------------------------------------------
# all Marginals
let
    # MODELS = [ME_BOUNDED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    MODELS = [ME_EXPECTED, ME_HOMO]

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf], lw = 2)
        gparams = (xaxis = nothing, yaxis = nothing, grid = false)
        for rxn in InLP.RXNS
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in MODELS
                plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg", gparams...)        
        plot_pol!(p, ME_EXPECTED, Vl, D, ϵ, τ; sparams...)
        plot_pol!(p, ME_HOMO, Vl, D, ϵ, τ; sparams...)
        plot_pol!(p, ME_BOUNDED, Vl, D, ϵ, τ; sparams...)
        push!(ps, p)
        
        mysavefig(ps, "All_Marginals"; Vl, D, ϵ)
    end
end

## ----------------------------------------------------------------------------
# vatp, vg marginals
let
    # MODELS = [ME_BOUNDED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    MODELS = [ME_BOUNDED, ME_EXPECTED, ME_CUTTED, ME_HOMO]

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf], lw = 3)
        gparams = (;grid = false)
        for rxn in ["gt", "vatp", "resp", "lt", "ldh"]
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in MODELS
                plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        p = plot(;title = "dynamic", xlabel = "", ylabel = "conc", gparams...)        
        M = idxdat([ME_EXPECTED, :M], Vl, D, ϵ, τ) 
        bar!(p, ["sg", "sl"], [M.sg, M.sl]; label = "")
        push!(ps, p)
        mysavefig(ps, "Marginals_v2"; Vl, D, ϵ)
    end
end

## ----------------------------------------------------------------------------
# Dev
let
    # Vl, D, ϵ, τ = EXP_PARAMS |> collect |> rand
    Vl, D, ϵ, τ = (0.0, 0.003, 0.01, 0.0)
    status = MINDEX[:STATUS, Vl, D, ϵ, τ]
    
    DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
    
    MODsym = ME_CUTTED
    M = idxdat([MODsym, :M], Vl, D, ϵ, τ) 
    MEM2s = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
    net = M.net
    net.ub[7] = 0.0
    L, U = InLP.fva(M.net)
    net.lb .= L; net.ub .= U

    δ = 0.08
    y = InLP.Y # atp/biomass yield
    maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)
    MEMs = InLP.get_marginals(maxentf(0.0), M; δ)
    
    for (rxni, rxn) in net.rxns |> enumerate
        DYav = InLP.av(DyMs[rxn])
        MEav = InLP.av(MEMs[rxn])
        ME2av = InLP.av(MEM2s[rxn])
        @info(rxn, (Vl, D, ϵ, τ), status, DYav, MEav, ME2av, net.lb[rxni], net.ub[rxni])
        println()
    end
end

## ----------------------------------------------------------------------------
# let
#     Vl = MINDEX[:Vls] |> first
#     τ = MINDEX[:τs] |> first
#     D = (MINDEX[:Ds] |> sort)[6]
#     ϵ = MINDEX[:ϵs][2]
#     MODsym = FBA_BOUNDED   
#     M = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
#     # @show MINDEX[:STATUS, Vl, D, ϵ, τ]
#     # vatp_range, vg_ranges = idxdat([ME_EXPECTED, :Ms], Vl, D, ϵ, τ)
#     # dat = deserialize(MINDEX[:DFILE, Vl, D, ϵ, τ]);
#     # dat[[ME_EXPECTED, :POL]...] |> keys |> collect
# end

## ----------------------------------------------------------------------------
# selected Marginals
let

    Vl = MINDEX[:Vls] |> first
    τ = MINDEX[:τs] |> first
    D = (MINDEX[:Ds] |> sort)[5]
    ϵs = MINDEX[:ϵs] |> sort
    sparams =(;alpha = 0.8, lw = 10, ylim = [0.0, Inf])
    gparams = (;grid = false)
    
    MODELS = [ME_BOUNDED, ME_EXPECTED, ME_CUTTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    ps = Plots.Plot[]
    for ϵ in ϵs
        for rxn in ["gt", "vatp"]
            MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
            @info("Doing", (rxn, Vl, D, ϵ, τ))
            p = plot(;title = "ϵ = $ϵ", 
                xlabel = rxn == "gt" ? "vg" : rxn, 
                ylabel = "prob", gparams...
            )
            
            @time begin
                plot_marginals!(p, MODELS, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end
    end
    layout = 4, 2
    mysavefig(ps, "dyn_vs_model_marginals"; layout, Vl, D, τ)
end

## ----------------------------------------------------------------------------
# beta vs eps
let
    ϵs = MINDEX[:ϵs] |> sort
    colors = Plots.distinguishable_colors(length(MINDEX[:Ds]))
    colors = Dict(D => c for (D, c) in zip(MINDEX[:Ds], colors))
    p = plot(;xlabel = "beta", ylabel = "ϵ")
    sparams = (;alpha = 0.5, ms = 6)
    exp_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    for (Vl, D, τ) in exp_params
        ϵ_ser = []
        beta_ser = []
        for ϵ in ϵs
            MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
            beta = idxdat([:ME_EXPECTED, :beta0], Vl, D, ϵ, τ)
            push!(ϵ_ser, ϵ)
            push!(beta_ser, beta)
        end
        scatter!(p, beta_ser, ϵ_ser; label = "", color = colors[D], sparams...)
    end
    mysavefig(p, "$(ME_EXPECTED)_beta_vs_eps_D_colored")
end

## ----------------------------------------------------------------------------
# D vs vatp
let
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort
    τ = MINDEX[:τs] |> first
    Vl = MINDEX[:Vls] |> first

    p = plot(;title = "dynamic stst", xlabel = "D", ylabel = "vatp")
    for ϵ in ϵs
        @info("Doing", ϵ)
        xs, ys = [], []
        for D in Ds
            MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
            DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
            
            vatp = InLP.av(DyMs["vatp"])
            push!(xs, D) 
            push!(ys, vatp) 
        end
        plot!(p, xs, ys; label = ϵ, alpha = 0.5, lw = 3)
    end
    mysavefig(p, "mu_vs_eps") 
end

## ----------------------------------------------------------------------------
# vatp corrs
let
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort
    τ = MINDEX[:τs] |> first
    Vl = MINDEX[:Vls] |> first

    MODELS = [ME_BOUNDED, ME_EXPECTED, ME_CUTTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]

    for flx in ["vatp", "gt"]
        ps = Dict(
            MODsym => plot(;title = string("dynamic stst: ", MODsym), 
                xlabel = "dyn $flx", ylabel = "model $flx"
            ) 
            for MODsym in MODELS
        )
        for ϵ in ϵs
            @info("Doing", ϵ)
            for D in Ds
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                dym_vatp = InLP.av(DyMs[flx])
                
                for MODsym in MODELS
                    Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                    m_vatp = InLP.av(Ms[flx])
                    color = MOD_COLORS[MODsym]
                    scatter!(ps[MODsym], [dym_vatp], [m_vatp]; color, 
                        label = "", alpha = 0.5, m = 8
                    )
                end
            end
        end
        ps = collect(values(ps))
        mysavefig(ps, "$(flx)_correlation") 
    end
end

## ----------------------------------------------------------------------------
# Steady State Model Dynamic correlation
let
    f(x) = log10(abs(x) + 1e-8)

    # PARAMS
    ϵs = MINDEX[:ϵs]
    sim_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    sim_params = collect(sim_params)[1:5:end]
    models = [ME_BOUNDED, ME_EXPECTED, ME_CUTTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    
    ps = Plots.Plot[]
    for ϵ in ϵs |> sort

        for  MODsym in models
            
            p = plot(;title = string(MODsym, " ϵ: ", ϵ), 
                xlabel = "dym flxs", ylabel = "model flxs", 
                legend = :topleft
            )
            
            color = MOD_COLORS[MODsym]
            arr = Vector{Float64}(undef, length(sim_params) * length(InLP.RXNS))
            DYN_flxs, DYN_errs = arr, copy(arr)
            M_flxs, M_errs = copy(arr), copy(arr)
            
            @info("Doing", ϵ, MODsym)
            
            for (i, (Vl, D, τ)) in sim_params |> enumerate
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                
                for rxn in InLP.RXNS
                    DYN_flx = InLP.av(DyMs[rxn])
                    DYN_err = InLP.va(DyMs[rxn]) |> sqrt
                    ME_flx = InLP.av(Ms[rxn])
                    ME_err = InLP.va(Ms[rxn]) |> sqrt
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
    rows = Int(round(length(ps) / length(models), RoundUp))
    cols = length(models)
    layout = rows, cols
    mysavefig(ps, "flxs_corr"; layout) 
end
