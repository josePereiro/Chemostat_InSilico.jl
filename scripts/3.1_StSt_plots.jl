import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic
        
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

const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_OPEN_G_BOUNDED       = :ME_Z_OPEN_G_BOUNDED        # 

const ME_Z_EXPECTED_G_OPEN      = :ME_Z_EXPECTED_G_OPEN       # Match ME and Dy biom average
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED   # 
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 

const ME_Z_FIXXED_G_OPEN        = :ME_Z_FIXXED_G_OPEN         # Fix biom around observed
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

const FBA_Z_OPEN_G_OPEN       = :FBA_Z_OPEN_G_OPEN
const FBA_Z_OPEN_G_BOUNDED    = :FBA_Z_OPEN_G_BOUNDED
const FBA_Z_FIXXED_G_OPEN     = :FBA_Z_FIXXED_G_OPEN 
const FBA_Z_FIXXED_G_BOUNDED  = :FBA_Z_FIXXED_G_BOUNDED

ALL_MODELS = [
    # ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
    # ME_Z_EXPECTED_G_OPEN, 
    ME_Z_EXPECTED_G_BOUNDED, 
    ME_Z_EXPECTED_G_MOVING,
    # ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_EXPECTED,
    # FBA_Z_OPEN_G_OPEN, FBA_Z_OPEN_G_BOUNDED, 
    # FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED
]

MOD_COLORS = Dict(
    ME_Z_OPEN_G_OPEN        => :brown, 
    ME_Z_OPEN_G_BOUNDED     => :orange, 
    ME_Z_EXPECTED_G_OPEN    => :red, 
    ME_Z_EXPECTED_G_BOUNDED => :blue,
    ME_Z_EXPECTED_G_MOVING => :purple,
    ME_Z_FIXXED_G_OPEN      => :green, 
    ME_Z_FIXXED_G_BOUNDED   => :pink, 
    ME_Z_EXPECTED_G_EXPECTED   => :violet, 

    FBA_Z_OPEN_G_OPEN       => :dot, 
    FBA_Z_OPEN_G_BOUNDED    => :dot, 
    FBA_Z_FIXXED_G_OPEN     => :dot, 
    FBA_Z_FIXXED_G_BOUNDED  => :dot,
)

MOD_LS = Dict(
    ME_Z_OPEN_G_OPEN        => :dash, 
    ME_Z_OPEN_G_BOUNDED     => :dash, 
    ME_Z_EXPECTED_G_OPEN    => :dash, 
    ME_Z_EXPECTED_G_BOUNDED => :dash,
    ME_Z_EXPECTED_G_MOVING => :dash,
    ME_Z_FIXXED_G_OPEN      => :dash, 
    ME_Z_FIXXED_G_BOUNDED   => :dash, 
    ME_Z_EXPECTED_G_EXPECTED   => :dash, 

    FBA_Z_OPEN_G_OPEN       => :dot, 
    FBA_Z_OPEN_G_BOUNDED    => :dot, 
    FBA_Z_FIXXED_G_OPEN     => :dot, 
    FBA_Z_FIXXED_G_BOUNDED  => :dot,
)
 

## ----------------------------------------------------------------------------
# MINDEX
# MDAT[MODsym, :M, Vl, D, ϵ, τ]
# MDAT[MODsym, :Ms, Vl, D, ϵ, τ]
# MDAT[MODsym, :beta_biom, Vl, D, ϵ, τ]
# MDAT[:STATUS, Vl, D, ϵ]
MINDEX_FILE = Dyn.procdir("marg_dat_index.bson")
MINDEX = UJL.load_data(MINDEX_FILE)
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...)
idxdat(dk, indexks...) = Dyn.idxdat(MINDEX, dk, indexks...)

## ----------------------------------------------------------------------------
# PLOTS
mysavefig(p, pname; params...) = 
    Dyn.mysavefig(p, pname, Dyn.plotsdir(), fileid; params...)

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

function plot_marginals!(p, MODsyms, rxn, Vl, D, ϵ, τ; draw_av = true, sparams...)

    DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
    plot!(p, DyMs[rxn]; label = "", sparams..., color = :black)
    if draw_av
        av = Dyn.av(DyMs[rxn])
        vline!(p, [av]; label = "", sparams..., 
            color = :black, ls = :solid, lw = 5, alpha = 0.6
        )
    end
    
    # Marginals
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        plot!(p, Ms[rxn]; label = "", color, ls, sparams...)
        if draw_av
            av = Dyn.av(Ms[rxn])
            vline!(p, [av]; label = "", color, sparams..., 
                ls = :solid, lw = 5, alpha = 0.6
            )
        end
    end
    return p
end
plot_marginals!(p, MODsyms::Symbol, rxn, Vl, D, ϵ, τ; sparams...) = 
    plot_marginals!(p, [MODsyms], rxn, Vl, D, ϵ, τ; sparams...)


## ----------------------------------------------------------------------------
# Check dyn quality
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
        av_z = Dyn.av(DyMs["biom"])
        μ = av_z - M.σ
        growth_bal = (μ - D)/D
        biom_prod = μ * M.X
        bioms_drain = D * M.X

        av_vg = Dyn.av(DyMs["gt"])
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

## ----------------------------------------------------------------------------
# Error
let
    Vl = MINDEX[:Vls] |> first
    τ = MINDEX[:τs] |> first
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort

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
                
                for rxn in Dyn.RXNS
                    DYN_flx = Dyn.av(DyMs[rxn])
                    ME_flx = Dyn.av(Ms[rxn])
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

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf], lw = 4)
        gparams = (;xaxis = nothing, yaxis = nothing, grid = false)
        for rxn in ["vatp", "gt", "biom"]
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in ALL_MODELS
                plot_marginals!(p, MODsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        # LEGEND
        p = plot(;title = "Legend", xaxis = nothing, yaxis = nothing)
        for  (i, MODsym) in ALL_MODELS |> enumerate
            ls = MOD_LS[MODsym] 
            color = MOD_COLORS[MODsym]
            plot!(p, fill(i, 3); label = string(MODsym), color, ls, lw = 8)        
        end
        push!(ps, p)
        
        mysavefig(ps, "All_Marginals"; Vl, D, ϵ)
    end
end

## ----------------------------------------------------------------------------
# let

#     params = EXP_PARAMS |> collect
#     δ = 0.08
#     LP_cache = nothing

#     @threads for (Vl, D, ϵ, τ) in params
#         status = MINDEX[:STATUS, Vl, D, ϵ, τ]
#         status != :stst && continue
        
#         # LOAD DATA
#         cfile = MINDEX[:DFILE, Vl, D, ϵ, τ]

#         @info("Doing", 
#             (Vl, D, ϵ, τ), δ
#         )

#         # RECOMPUTE MARGINALS
#         M0 = idxdat([:M0], Vl, D, ϵ, τ)
#         LP_cache = isnothing(LP_cache) ? Dyn.vgvatp_cache(M0) : LP_cache
#         MDAT = deserialize(cfile)

#         @info("Dynamic") 
#         MDAT[:DyMs] = Dyn.get_marginals(M0; δ, LP_cache, verbose = false) do vatp, vg
#             M0.Xb[vatp][vg] / M0.X
#         end

#         # EP
#         for MODsym in [
#                     ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
#                     ME_Z_EXPECTED_G_OPEN, ME_Z_EXPECTED_G_BOUNDED,
#                     ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
#                     ME_Z_EXPECTED_G_EXPECTED
#                 ]
#             M = idxdat([MODsym, :M], Vl, D, ϵ, τ)
#             beta_biom = idxdat([MODsym, :beta_biom], Vl, D, ϵ, τ)
#             beta_vg = idxdat([MODsym, :beta_vg], Vl, D, ϵ, τ)
            
#             @info("EP", MODsym, beta_biom, beta_vg)
#             z(vatp, vg) = LP_cache[vatp][vg][M.obj_idx]
#             MDAT[MODsym, :Ms] = Dyn.get_marginals(M; δ, LP_cache, verbose = false) do vatp, vg
#                 exp((beta_biom * z(vatp, vg)) + (beta_vg * vg))
#             end
#         end

#         # FBA
#         for MODsym in [
#                     FBA_Z_OPEN_G_OPEN, FBA_Z_OPEN_G_BOUNDED, 
#                     FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED
#                 ]
#             M = idxdat([MODsym, :M], Vl, D, ϵ, τ)
            
#             @info("FBA", MODsym)
#             # FBA
#             # Find maximum feasible vatp
#             vatp_range, vg_ranges = Dyn.vatpvg_ranges(M)
#             max_vatp = -Inf
#             min_vg = Inf
#             # (max_vatp, min_vg) will maximize the yield
#             for (vatpi, vatp) in vatp_range |> enumerate
#                 vg_range = vg_ranges[vatpi]
#                 isempty(vg_range) && continue
#                 if vatp > max_vatp
#                     max_vatp = vatp
#                     min_vg = minimum(vg_range)
#                 end
#             end
#             @assert !isinf(max_vatp)

#             fbaf(vatp, vg) = (vatp == max_vatp && vg == min_vg) ? 1.0 : 0.0
#             MDAT[MODsym, :Ms] = Dyn.get_marginals(fbaf, M; δ, LP_cache, verbose = false)
#         end


#         # SAVING
#         serialize(cfile, MDAT)
#         Dyn.idxdat(;emptytcache = true)
        
#         # # TEST
#         # rxn = "biom"
#         # Dyn.idxdat(;emptytcache = true)
#         # p = plot(;title = rxn, xlabel = "flx", ylabel = "prob")
#         # plot_marginals!(p, ALL_MODELS, rxn, Vl, D, ϵ, τ)
#         # mysavefig(p, "test")


#     end # for (Vl, D, ϵ, τ) in params
# end

## ----------------------------------------------------------------------------
# vatp, vg marginals v2
let
    MODELS = [
        # ME_Z_OPEN_G_OPEN, 
        # ME_Z_OPEN_G_BOUNDED, 
        # ME_Z_EXPECTED_G_OPEN, 
        ME_Z_EXPECTED_G_BOUNDED,
        :ME_Z_EXPECTED_G_MOVING,
        ME_Z_EXPECTED_G_EXPECTED,
        # ME_Z_FIXXED_G_OPEN, 
        # ME_Z_FIXXED_G_BOUNDED, 
        # FBA_Z_OPEN_G_OPEN, 
        # FBA_Z_OPEN_G_BOUNDED, 
        # FBA_Z_FIXXED_G_OPEN, 
        # FBA_Z_FIXXED_G_BOUNDED
    ]
    MODELS = ALL_MODELS

    # LEGEND PLOT
    leg_p = plot(;title = "Legend", xaxis = nothing, yaxis = nothing)
    for  (i, MODsym) in MODELS |> enumerate
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        plot!(leg_p, fill(length(MODELS) - i, 3); label = string(MODsym), color, ls, lw = 8)        
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
            for  MODsym in MODELS
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
#             beta = idxdat([:ME_EXPECTED, :beta_biom], Vl, D, ϵ, τ)
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
            
            vatp = Dyn.av(DyMs["vatp"])
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
# \BETA vs D vs \EPS
let
    MODsym = ME_Z_EXPECTED_G_BOUNDED

    dat_pool = Dict()
    for (Vl, D, ϵ, τ) in EXP_PARAMS
        dat = get!(dat_pool, ϵ, Dict())
        get!(dat, :Ds, [])
        get!(dat, :betas, [])

        status = MINDEX[:STATUS, Vl, D, ϵ, τ] 
        
        @info("Doing", (Vl, D, ϵ, τ), MODsym, status)
        status != :stst && continue
        beta = idxdat([MODsym, :beta_biom], Vl, D, ϵ, τ)
        # beta = rand() # Test

        push!(dat[:Ds], D)
        push!(dat[:betas], beta)
    end

    # PLOTTING D vs BETA
    p = plot(;xlabel = "D", ylabel = "beta")
    for (eps, dat) in dat_pool
        color = Gray(eps * 0.8)
        scatter!(p, dat[:Ds], dat[:betas]; 
            label = "", color, m = 8, alpha = 0.7
        )
    end
    mysavefig(p, "beta_vs_D"; MODsym) 

    # PLOTTING EPS vs BETA
    p = plot(;xlabel = "eps", ylabel = "beta")
    for (eps, dat) in dat_pool
        epss = fill(eps, length(dat[:betas]))
        scatter!(p, epss, dat[:betas]; 
            label = "", color = :black, m = 8, alpha = 0.7
        )
    end
    mysavefig(p, "beta_vs_eps"; MODsym) 
end


## ----------------------------------------------------------------------------
let
    # LEGEND
    leg_p = plot(;title = "legend", xaxis = nothing, ylabel = "ϵ")
    for ϵ in MINDEX[:ϵs]
        color = Gray(ϵ * 0.8)
        plot!(leg_p, fill(ϵ, 10); color, label = string("ϵ: ", ϵ), lw = 8)
    end
    mysavefig(leg_p, "eps_legend") 
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