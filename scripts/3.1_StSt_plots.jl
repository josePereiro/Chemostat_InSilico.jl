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

const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_OPEN_G_BOUNDED     = :ME_Z_OPEN_G_BOUNDED        # 
const ME_Z_EXPECTED_G_OPEN    = :ME_Z_EXPECTED_G_OPEN       # Match ME and Dy biom average
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_FIXXED_G_OPEN      = :ME_Z_FIXXED_G_OPEN         # Fix biom around observed
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

const FBA_Z_OPEN_G_OPEN       = :FBA_Z_OPEN_G_OPEN
const FBA_Z_OPEN_G_BOUNDED    = :FBA_Z_OPEN_G_BOUNDED
const FBA_Z_FIXXED_G_OPEN     = :FBA_Z_FIXXED_G_OPEN 
const FBA_Z_FIXXED_G_BOUNDED  = :FBA_Z_FIXXED_G_BOUNDED

ALL_MODELS = [
    ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
    ME_Z_EXPECTED_G_OPEN, ME_Z_EXPECTED_G_BOUNDED,
    ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
    FBA_Z_OPEN_G_OPEN, FBA_Z_OPEN_G_BOUNDED, 
    FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED
]

MOD_COLORS = let
    colors = Plots.distinguishable_colors(length(ALL_MODELS))
    Dict(mod => color for (mod, color) in zip(ALL_MODELS, colors))
end

MOD_LS = Dict(
    ME_Z_OPEN_G_OPEN        => :dash, 
    ME_Z_OPEN_G_BOUNDED     => :dash, 
    ME_Z_EXPECTED_G_OPEN    => :dash, 
    ME_Z_EXPECTED_G_BOUNDED => :dash,
    ME_Z_FIXXED_G_OPEN      => :dash, 
    ME_Z_FIXXED_G_BOUNDED   => :dash, 

    FBA_Z_OPEN_G_OPEN       => :dot, 
    FBA_Z_OPEN_G_BOUNDED    => :dot, 
    FBA_Z_FIXXED_G_OPEN     => :dot, 
    FBA_Z_FIXXED_G_BOUNDED  => :dot,
)
 

## ----------------------------------------------------------------------------
# MINDEX
# MDAT[MODsym, :M, Vl, D, ϵ, τ]
# MDAT[MODsym, :Ms, Vl, D, ϵ, τ]
# MDAT[MODsym, :beta0, Vl, D, ϵ, τ]
# MDAT[:STATUS, Vl, D, ϵ]
MINDEX = UJL.load_data(InCh.MARGINALS_INDEX_FILE)
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...)
idxdat(dk, indexks...) = InLP.idxdat(MINDEX, dk, indexks...)

# ## ----------------------------------------------------------------------------
# let
#     for (Vl, D, ϵ, τ) in collect(EXP_PARAMS)
#         status = MINDEX[:STATUS, Vl, D, ϵ, τ]
#         dfile = MINDEX[:DFILE, Vl, D, ϵ, τ]
#         absfile = joinpath(InCh.PROJECT_DIR, dfile)
#         @info("Doing", 
#             (Vl, D, ϵ, τ),
#             status, dfile,
#             isfile(absfile)
#         )
#         @assert status != :stst || (status == :stst && isfile(absfile))
#     end
# end

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


## ----------------------------------------------------------------------------
# Check dyn constraints
let

    biom_prods, bioms_drains = [], []
    glc_ins, glc_ups, glc_drains = [], [], []
    # LP_cache = nothing
    for (Vl, D, ϵ, τ) in collect(EXP_PARAMS)
        
        # LOAD
        status = MINDEX[:STATUS, Vl, D, ϵ, τ]
        @info("Doing", (Vl, D, ϵ, τ), status)
        status != :stst && continue
        println()

        DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        M = idxdat([:M0], Vl, D, ϵ, τ)
        
        # COMPUTE
        f(vatp, vg) = M.Xb[vatp][vg] / M.X
        av_z = InLP.av(DyMs["biom"])
        μ = av_z - M.σ
        growth_bal = (μ - D)/D
        biom_prod = μ * M.X
        bioms_drain = D * M.X

        av_vg = InLP.av(DyMs["gt"])
        cD_X = M.cg * M.D / M.X
        glc_bal = (-av_vg * M.X + (M.cg - M.sg) * M.D) / (M.cg * M.D)
        glc_in = M.cg * M.D
        glc_drain = M.sg * M.D
        glc_up = av_vg * M.X
        
        # PUSH
        @info("Growth", av_z, μ, M.σ, D, growth_bal)
        @info("GLC uptake", av_vg, cD_X, M.sg, glc_bal)
        println()

        push!(biom_prods, biom_prod)
        push!(bioms_drains, bioms_drain)
        push!(glc_ins, glc_in)
        push!(glc_drains, glc_drain)
        push!(glc_ups, glc_up)
    end

    # PLOT
    biom_bal_p = plot(;title = "Biomass balance", 
        xlabel = "simulation", ylabel = "rate"
    )
    plot!(biom_bal_p, biom_prods; label = "prods", lw = 3)
    plot!(biom_bal_p, bioms_drains; label = "drains", lw = 3)

    biom_corr_p = plot(title = "Biomass balance correlation", 
       xlabel = "production", ylabel = "drain"
    )
    scatter!(biom_corr_p, biom_prods, bioms_drains; label = "", m = 8)
    vals = sort([biom_prods; bioms_drains])
    plot!(biom_corr_p, vals, vals; label = "", color = :black, ls = :dash, alpha = 0.7)
    
    glc_bal_p = plot(;title = "Glc balance", 
        xlabel = "simulation", ylabel = "rate"
    )
    plot!(glc_bal_p, glc_ins; label = "input", lw = 3)
    plot!(glc_bal_p, glc_ups; label = "uptake", lw = 3)
    plot!(glc_bal_p, glc_drains; label = "drain", lw = 3)
    plot!(glc_bal_p, glc_drains .+ glc_ups; label = "up + drain", lw = 3)

    glc_corr_b = plot(title = "Glc balance correlation", 
        xlabel = "production", ylabel = "up + drain"
    )
    scatter!(glc_corr_b,  glc_ins, glc_drains .+ glc_ups; label = "", m = 8)
    vals = sort([glc_ins; glc_drains .+ glc_ups])
    plot!(glc_corr_b, vals, vals; label = "", color = :black, ls = :dash, alpha = 0.7)

    ps = Plots.Plot[biom_bal_p, biom_corr_p, glc_bal_p, glc_corr_b]
    mysavefig(ps, "dyn_balances")

end

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

    # MODELS = [ME_BOUNDED, ME_CUTTED, ME_EXPECTED, ME_Z_OPEN_G_OPEN, FBA_BOUNDED, FBA_OPEN]
    ps = Plots.Plot[]
    M = -Inf
    for MODsym in ALL_MODELS
        p = plot(;tile = string("Error", MODsym), 
            xlabel = "log ϵ", ylabel = "maximum err"
        )
        xs, ys, yerrs = [], [], []
        for ϵ in ϵs
            errs = []
            @info("Doing", ϵ, MODsym); println()
            for D in Ds
                status = MINDEX[:STATUS, Vl, D, ϵ, τ]
                status != :stst && continue
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
            max_err = maximum(errs)
            push!(ys, max_err)
            M = max(M, max_err)

        end # for ϵ in ϵs

        noise = xs .* 0.1 .* rand.()
        params = (;alpha = 0.8, color = MOD_COLORS[MODsym])
        plot!(p, log10.(xs .+ noise), ys; 
            label = "", lw = 3, ls = :dash, params...
        )
        scatter!(p, log10.(xs .+ noise), ys; # yerr = yerrs, 
            ms = 8, params..., label = string(MODsym), legend = :topleft
        )
        push!(ps, p)
    end #  for MODsym 

    mysavefig(ps, "eps_vs_err")
    
    for p in ps
        plot!(p; ylim = [0.0, M * 1.1])
    end
    mysavefig(ps, "eps_vs_err_equal_scale")
    

end

## ----------------------------------------------------------------------------
# all Marginals
let
    # MODELS = [ME_BOUNDED, ME_EXPECTED, ME_Z_OPEN_G_OPEN, FBA_BOUNDED, FBA_OPEN]
    # MODELS = [ME_EXPECTED, ME_CUTTED, ME_Z_OPEN_G_OPEN, FBA_BOUNDED]

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf], lw = 4)
        gparams = (xaxis = nothing, yaxis = nothing, grid = false)
        for rxn in ["vatp", "gt", "biom"]
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in ALL_MODELS
                plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        # LEGEND
        p = plot(;title = "Legend", xaxis = nothing, yaxis = nothing)
        for  MODsym in ALL_MODELS
            ls = MOD_LS[MODsym] 
            color = MOD_COLORS[MODsym]
            plot!(p, rand(3); label = string(MODsym), color, ls, lw = 8)        
        end
        push!(ps, p)
        
        mysavefig(ps, "All_Marginals"; Vl, D, ϵ)
    end
end

## ----------------------------------------------------------------------------
# vatp, vg marginals
let
    MODELS = [
        ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
        ME_Z_EXPECTED_G_OPEN, ME_Z_EXPECTED_G_BOUNDED,
        # # ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
        # FBA_Z_OPEN_G_OPEN, FBA_Z_OPEN_G_BOUNDED, 
        # FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED
    ]

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf], lw = 4)
        gparams = (;grid = false)
        for rxn in ["gt", "biom", "resp", "lt"]
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in MODELS
                plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        M = idxdat([:M0], Vl, D, ϵ, τ) 
        
        # POLYTOPE
        p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg")
        InLP.plot_polborder!(p, M)
        InLP.plot_poldist!(p, M)
        DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        vatp_av = InLP.av(DyMs["vatp"]) 
        vg_av = InLP.av(DyMs["gt"]) 

        kwargs = (;alpha = 0.5, label = "", color = :black)
        hline!(p, [vg_av]; lw = 5, ls = :dash, kwargs...)
        hline!(p, [M.cg * M.D / M.X]; lw = 5, ls = :dot, kwargs...)
        vline!(p, [vatp_av]; lw = 5, ls = :dash, kwargs...)
        push!(ps, p)

        # LEGEND
        p = plot(;title = "Legend", xaxis = nothing, yaxis = nothing)
        for  MODsym in MODELS
            ls = MOD_LS[MODsym] 
            color = MOD_COLORS[MODsym]
            plot!(p, rand(3); label = string(MODsym), color, ls, lw = 8)        
        end
        plot!(p, rand(3); label = "DYNAMIC", color = :black, lw = 8)        
        push!(ps, p)

        mysavefig(ps, "Marginals_v2"; Vl, D, ϵ, M.sg, M.sl)
    end
end

## ----------------------------------------------------------------------------
# # Dev
# let
#     # Vl, D, ϵ, τ = EXP_PARAMS |> collect |> rand
#     Vl, D, ϵ, τ = (0.0, 0.003, 0.01, 0.0)
#     status = MINDEX[:STATUS, Vl, D, ϵ, τ]
    
#     DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
    
#     MODsym = ME_CUTTED
#     M = idxdat([MODsym, :M], Vl, D, ϵ, τ) 
#     MEM2s = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
#     net = M.net
#     net.ub[7] = 0.0
#     L, U = InLP.fva(M.net)
#     net.lb .= L; net.ub .= U

#     δ = 0.08
#     y = InLP.Y # atp/biomass yield
#     maxentf(beta) = (vatp, vg) -> exp(beta * vatp/y)
#     MEMs = InLP.get_marginals(maxentf(0.0), M; δ)
    
#     for (rxni, rxn) in net.rxns |> enumerate
#         DYav = InLP.av(DyMs[rxn])
#         MEav = InLP.av(MEMs[rxn])
#         ME2av = InLP.av(MEM2s[rxn])
#         @info(rxn, (Vl, D, ϵ, τ), status, DYav, MEav, ME2av, net.lb[rxni], net.ub[rxni])
#         println()
#     end
# end

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

## ----------------------------------------------------------------------------
# # beta vs eps
# let
#     ϵs = MINDEX[:ϵs] |> sort
#     colors = Plots.distinguishable_colors(length(MINDEX[:Ds]))
#     colors = Dict(D => c for (D, c) in zip(MINDEX[:Ds], colors))
#     p = plot(;xlabel = "beta", ylabel = "ϵ")
#     sparams = (;alpha = 0.5, ms = 6)
#     exp_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
#     for (Vl, D, τ) in exp_params
#         ϵ_ser = []
#         beta_ser = []
#         for ϵ in ϵs
#             MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
#             beta = idxdat([:ME_EXPECTED, :beta0], Vl, D, ϵ, τ)
#             push!(ϵ_ser, ϵ)
#             push!(beta_ser, beta)
#         end
#         scatter!(p, beta_ser, ϵ_ser; label = "", color = colors[D], sparams...)
#     end
#     mysavefig(p, "$(ME_EXPECTED)_beta_vs_eps_D_colored")
# end

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
# Steady State Model Dynamic correlation
let
    f(x) = log10(abs(x) + 1e-8)

    # PARAMS
    ϵs = MINDEX[:ϵs]
    sim_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    sim_params = collect(sim_params)[1:5:end]
    # models = [ME_BOUNDED, ME_EXPECTED, ME_CUTTED, ME_Z_OPEN_G_OPEN, FBA_BOUNDED, FBA_OPEN]
    
    ps = Plots.Plot[]
    for ϵ in ϵs |> sort

        for  MODsym in ALL_MODELS
            
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
    rows = Int(round(length(ps) / length(ALL_MODELS), RoundUp))
    cols = length(ALL_MODELS)
    layout = rows, cols
    mysavefig(ps, "flxs_corr"; layout) 
end


## ----------------------------------------------------------------------------
# flx corrs
let
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort
    τ = MINDEX[:τs] |> first
    Vl = MINDEX[:Vls] |> first

    # gps = Dict()
    color_pool = Dict(
        ϵ => c for (ϵ, c) in 
        zip(ϵs, Plots.distinguishable_colors(length(ϵs)))
    )
    for flx in ["vatp", "gt"]           
        colors, ys, xs = Dict(), Dict(), Dict()
        for ϵ in ϵs
            @info("Doing", ϵ)
            for D in Ds
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                dym_vatp = InLP.av(DyMs[flx])
                # dym_vatp = rand() # Test

                for MODsym in ALL_MODELS
                    get!(xs, MODsym, [])
                    get!(ys, MODsym, [])
                    get!(colors, MODsym, [])
                    Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
                    m_vatp = InLP.av(Ms[flx])
                    # m_vatp = rand() # Test

                    push!(xs[MODsym], dym_vatp)
                    push!(ys[MODsym], m_vatp)
                    push!(colors[MODsym], color_pool[ϵ])
                end
            end
        end

        # @show color
        lps = Plots.Plot[]
        for MODsym in ALL_MODELS
            lp = plot(; title = string("dynamic stst: ", MODsym), 
                xlabel = "dyn $flx", ylabel = "model $flx"
            ) 

            color = colors[MODsym]
            scatter!(lp, xs[MODsym], ys[MODsym]; color, 
                label = "", m = 8
            )
            vals = sort([xs[MODsym]; ys[MODsym]])
            plot!(lp, vals, vals; label = "", color = :black, 
                alpha = 0.8, lw = 3, ls = :dash
            )
            push!(lps,lp)
            
        end

        # legend
        p = plot(;title = "legend", xaxis = nothing, ylabel = "ϵ")
        for ϵ in ϵs
            color = color_pool[ϵ]
            plot!(p, fill(ϵ, 10); color, label = string("ϵ: ", ϵ), lw = 8)
        end
        push!(lps,p)

        mysavefig(lps, "$(flx)_correlation") 
    end
end