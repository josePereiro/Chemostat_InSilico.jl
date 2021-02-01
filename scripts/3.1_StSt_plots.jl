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
const STST_POL = :STST_POL         # Use polytope computed from chemostat stst assertion
const DYN_POL = :DYN_POL           # Use dynamic polytope

const ME_HOMO = :ME_HOMO           # Do not use extra constraints
const ME_EXPECTED = :ME_EXPECTED   # Match ME and Dy biom average
const ME_BOUNDED = :ME_BOUNDED     # Fix biom around observed

const FBA_OPEN = :FBA_OPEN
const FBA_BOUNDED = :FBA_BOUNDED

MOD_COLORS = Dict(
    ME_HOMO => :red, 
    ME_BOUNDED => :orange,
    ME_EXPECTED => :blue,
    FBA_BOUNDED => :green,
    FBA_OPEN => :purple,
)

POL_STYLE = Dict(
    STST_POL => :dot,
    DYN_POL => :dash,
)

## ----------------------------------------------------------------------------
# Marginals
# MDAT[MODsym, POLTsym, :M, Vl, D, ϵ, τ]
# MDAT[MODsym, POLTsym, :Ms, Vl, D, ϵ, τ]
# MDAT[MODsym, POLTsym, :beta0, Vl, D, ϵ, τ]
# MDAT[:STATUS, Vl, D, ϵ]
MINDEX = UJL.load_data(InCh.MARGINALS_INDEX_FILE)
POLTsym = STST_POL
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...)
idxdat(dk, indexks...) = InLP.idxdat(MINDEX, dk, indexks...)

## ----------------------------------------------------------------------------
# PLOTS
mysavefig(p, pname; params...) = 
    InLP.mysavefig(p, pname, InLP.DYN_FIGURES_DIR, fileid; params...)

## ----------------------------------------------------------------------------
# Plot functions
function plot_pol!(p, POLTsym, MODsym, Vl, D, ϵ, τ; sparams...)
    
    vatp_range, vg_ranges = idxdat([MODsym, POLTsym, :POL], Vl, D, ϵ, τ)
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
    params = (;label = "", ls = POL_STYLE[POLTsym], 
        alpha = 0.4, color, lw = 8, sparams...)
    plot!(p, [vatps], [vgLs]; params...)
    plot!(p, [vatps], [vgUs]; params...)
    p
end

function plot_marginals!(p, MODsyms, POLTsym, rxn, Vl, D, ϵ, τ; sparams...)

    ls = POL_STYLE[POLTsym]
    DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
    plot!(p, DyMs[rxn]; label = "", sparams..., color = :black)
    
    # Marginals
    for MODsym in MODsyms
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, POLTsym, :Ms], Vl, D, ϵ, τ)
        plot!(p, Ms[rxn]; label = "", color, sparams...)
    end
    return p
end
plot_marginals!(p, MODsyms::Symbol, POLTsym, rxn, Vl, D, ϵ, τ; sparams...) = 
    plot_marginals!(p, [MODsyms], POLTsym, rxn, Vl, D, ϵ, τ; sparams...)


## ----------------------------------------------------------------------------
# Error
let
    Vl = MINDEX[:Vls] |> first
    τ = MINDEX[:τs] |> first
    Ds = MINDEX[:Ds] |> sort
    ϵs = MINDEX[:ϵs] |> sort

    MODELS = [ME_BOUNDED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    POLsym = STST_POL

    p = plot(;tile = "Prediction error", xlabel = "log ϵ", ylabel = "err")
    for MODsym in MODELS
        xs, ys, yerrs = [], [], []
        for ϵ in ϵs

            errs = []
            @info("Doing", ϵ, MODsym); println()
            for D in Ds
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                Ms = idxdat([MODsym, POLsym, :Ms], Vl, D, ϵ, τ)
                
                for rxn in InLP.RXNS
                    DYN_flx = InLP.av(DyMs[rxn])
                    # DYN_err = InLP.va(DyMs[rxn]) |> sqrt
                    ME_flx = InLP.av(Ms[rxn])
                    # ME_err = InLP.va(Ms[rxn]) |> sqrt
                    (isnan(DYN_flx) || isnan(ME_flx)) && continue

                    err = ((DYN_flx - ME_flx)^2)/abs(DYN_flx)
                    push!(errs, err)
                end
            end

            push!(xs, ϵ)
            push!(ys, mean(errs))
            push!(yerrs, std(errs))

        end # for ϵ in ϵs

        noise = xs .* 0.1 .* rand.()
        scatter!(p, log10.(xs .+ noise), ys; yerr = yerrs, 
            alpha = 0.5, color = MOD_COLORS[MODsym], ms = 8,
            label = string(MODsym), legend = :topleft
        )
    end #  for MODsym 

    mysavefig(p, "eps_vs_err"; POLsym)

end

## ----------------------------------------------------------------------------
# all Marginals
let
    MODELS = [ME_BOUNDED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    return

    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
        ps = Plots.Plot[]
        sparams =(;ylim = [0.0, Inf])
        gparams = (xaxis = nothing, yaxis = nothing, grid = false)
        for rxn in InLP.RXNS
            p = plot(;title = rxn, xlabel = "flx", ylabel = "prob", gparams...)
            for  MODsym in MODELS
                plot_marginals!(p, MODsym, POLTsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end

        p = plot(;title = "polytope", xlabel = "vatp", ylabel = "vg", gparams...)        
        plot_pol!(p, POLTsym, ME_EXPECTED, Vl, D, ϵ, τ; sparams...)
        plot_pol!(p, POLTsym, ME_HOMO, Vl, D, ϵ, τ; sparams...)
        plot_pol!(p, POLTsym, ME_BOUNDED, Vl, D, ϵ, τ; sparams...)
        push!(ps, p)
        
        mysavefig(ps, "Marginals_$(POLTsym)_"; Vl, D, ϵ)
    end
end

## ----------------------------------------------------------------------------
# let
#     Vl = MINDEX[:Vls] |> first
#     τ = MINDEX[:τs] |> first
#     D = (MINDEX[:Ds] |> sort)[6]
#     ϵ = MINDEX[:ϵs][2]
#     MODsym = FBA_BOUNDED   
#     M = idxdat([MODsym, DYN_POL, :Ms], Vl, D, ϵ, τ)
#     # @show MINDEX[:STATUS, Vl, D, ϵ, τ]
#     # vatp_range, vg_ranges = idxdat([ME_EXPECTED, :STST_POL, :Ms], Vl, D, ϵ, τ)
#     # dat = deserialize(MINDEX[:DFILE, Vl, D, ϵ, τ]);
#     # dat[[ME_EXPECTED, :STST_POL, :POL]...] |> keys |> collect
# end

## ----------------------------------------------------------------------------
# selected Marginals
let

    Vl = MINDEX[:Vls] |> first
    τ = MINDEX[:τs] |> first
    D = (MINDEX[:Ds] |> sort)[9]
    ϵs = MINDEX[:ϵs] |> sort
    sparams =(;alpha = 0.5, lw = 10, ylim = [0.0, Inf])
    gparams = (;grid = false)
    
    MODELS = [ME_BOUNDED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]

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
                plot_marginals!(p, MODELS, POLTsym, rxn, Vl, D, ϵ, τ; sparams...)
            end
            push!(ps, p)
        end
    end
    layout = 4, 2
    mysavefig(ps, "dyn_vs_model_marginals"; layout, Vl, D, τ, POLTsym)
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
            beta = idxdat([:ME_EXPECTED, POLTsym, :beta0], Vl, D, ϵ, τ)
            push!(ϵ_ser, ϵ)
            push!(beta_ser, beta)
        end
        scatter!(p, beta_ser, ϵ_ser; label = "", color = colors[D], sparams...)
    end
    mysavefig(p, "$(ME_EXPECTED)_beta_vs_eps_D_colored")
end

## ----------------------------------------------------------------------------
# let
#     for (Vl, D, ϵ, τ) = EXP_PARAMS
#         MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
#         dat = idxdat([:DyMs], Vl, D, ϵ, τ)
#         @show length(dat)
#         break
#     end
# end

## ----------------------------------------------------------------------------
# Bound Correlation
let
    f(x) = x
    p = plot(;title = "Exchs Bounds Correlation", 
        xlabel = "dym bound", ylabel = "stst bound")
    l, u = Inf, -Inf

    STST_vgubs, DYN_vgubs = [], []
    STST_vlubs, DYN_vlubs = [], []
    for (Vl, D, ϵ, τ) in EXP_PARAMS
        MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue

        for  MODsym in [ME_BOUNDED, ME_EXPECTED, ME_HOMO]
            M = idxdat([MODsym, STST_POL, :M], Vl, D, ϵ, τ)
            push!(STST_vgubs, M.net.ub[M.vg_idx])
            push!(STST_vlubs, M.net.ub[M.vl_idx])

            M = idxdat([MODsym, DYN_POL, :M], Vl, D, ϵ, τ)
            push!(DYN_vgubs, M.net.ub[M.vg_idx])
            push!(DYN_vlubs, M.net.ub[M.vl_idx])
        end
    end

    xs = [DYN_vgubs; DYN_vlubs]
    ys = [STST_vgubs; STST_vlubs]
    l = minimum([l; xs; ys])            
    u = maximum([u; xs; ys])            
    plot!(p, f.([l, u]), f.([l, u]); label = "", ls = :dash, alpha = 0.8)
    scatter!(p, f.(xs), f.(ys); alpha = 0.5, color = :black, label = "")

    # saving
    mysavefig(p, "exchs_bounds_corr")

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

# ## ----------------------------------------------------------------------------
# # eps vs vg
# let
#     Ds = MINDEX[:Ds] |> sort
#     ϵs = MINDEX[:ϵs] |> sort
#     τ = MINDEX[:τs] |> first
#     Vl = MINDEX[:Vls] |> first

#     p = plot(;title = "dynamic stst", xlabel = "ϵ", ylabel = "vg")
#     for D in Ds
#         @info("Doing", D)
#         xs, ys = [], []
#         for ϵ in ϵs
#             MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
#             DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
            
#             av = InLP.av(DyMs["gt"])
#             push!(xs, ϵ) 
#             push!(ys, av) 
#         end
#         plot!(p, xs, ys; label = "", alpha = 0.5, lw = 3)
#     end
#     mysavefig(p, "eps_vs_vg") 
# end


## ----------------------------------------------------------------------------
# Steady State Model Dynamic correlation
let
    f(x) = log10(abs(x) + 1e-8)

    # PARAMS
    ϵs = MINDEX[:ϵs]
    sim_params = Iterators.product(MINDEX[[:Vls, :Ds, :τs]]...)
    sim_params = collect(sim_params)[1:5:end]
    models = [ME_BOUNDED, ME_EXPECTED, ME_HOMO, FBA_BOUNDED, FBA_OPEN]
    
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
            
            @info("Doing", ϵ, MODsym, STST_POL)
            
            for (i, (Vl, D, τ)) in sim_params |> enumerate
                MINDEX[:STATUS, Vl, D, ϵ, τ] != :stst && continue
                DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
                Ms = idxdat([MODsym, STST_POL, :Ms], Vl, D, ϵ, τ)
                
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
