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
    using ColorSchemes

    UJL.set_cache_dir(Dyn.cachedir())

end

# ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
minmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

# ----------------------------------------------------------------------------
# MINDEX
# MDAT[MODsym, :M, Vl, D, ϵ, τ]
# MDAT[MODsym, :Ms, Vl, D, ϵ, τ]
# MDAT[MODsym, :beta_biom, Vl, D, ϵ, τ]
# MDAT[:STATUS, Vl, D, ϵ]
MINDEX_FILE = Dyn.procdir("marg_dat_index.bson")
MINDEX = UJL.load_data(MINDEX_FILE)
EXP_PARAMS = Iterators.product(MINDEX[[:Vls, :Ds, :ϵs, :τs]]...)
idxdat(dk, indexks...) = Dyn.idxdat(MINDEX, dk, indexks...)

# ----------------------------------------------------------------------------

const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_OPEN_G_BOUNDED       = :ME_Z_OPEN_G_BOUNDED        # 

const ME_Z_EXPECTED_G_OPEN      = :ME_Z_EXPECTED_G_OPEN       # Match ME and Dy biom average
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED   # 
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_FULL_POLYTOPE          = :ME_FULL_POLYTOPE           # 

const ME_Z_FIXXED_G_OPEN        = :ME_Z_FIXXED_G_OPEN         # Fix biom around observed
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

const FBA_Z_OPEN_G_OPEN       = :FBA_Z_OPEN_G_OPEN
const FBA_Z_OPEN_G_BOUNDED    = :FBA_Z_OPEN_G_BOUNDED
const FBA_Z_FIXXED_G_OPEN     = :FBA_Z_FIXXED_G_OPEN 
const FBA_Z_FIXXED_G_BOUNDED  = :FBA_Z_FIXXED_G_BOUNDED

ME_MODELS = [
    ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
    ME_Z_EXPECTED_G_OPEN, 
    ME_Z_EXPECTED_G_BOUNDED, 
    ME_FULL_POLYTOPE,
    ME_Z_EXPECTED_G_MOVING,
    ME_Z_FIXXED_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
    ME_Z_EXPECTED_G_EXPECTED
]

FBA_MODELS = [
    FBA_Z_OPEN_G_OPEN, 
    FBA_Z_OPEN_G_BOUNDED, 
    FBA_Z_FIXXED_G_OPEN, 
    FBA_Z_FIXXED_G_BOUNDED
]

ALL_MODELS = [ME_MODELS; FBA_MODELS]

MOD_COLORS = Dict(
    :dyn => :black,
    ME_Z_OPEN_G_OPEN        => :brown, 
    ME_Z_OPEN_G_BOUNDED     => :orange, 
    ME_Z_EXPECTED_G_OPEN    => :red, 
    ME_Z_EXPECTED_G_BOUNDED => :blue,
    ME_FULL_POLYTOPE         => :red,
    ME_Z_EXPECTED_G_MOVING => :purple,
    ME_Z_FIXXED_G_OPEN      => :green, 
    ME_Z_FIXXED_G_BOUNDED   => :pink, 
    ME_Z_EXPECTED_G_EXPECTED   => :violet, 

    (FBA_MODELS .=> :blue)...
    # FBA_Z_OPEN_G_OPEN       => :blue, 
    # FBA_Z_OPEN_G_BOUNDED    => :blue, 
    # FBA_Z_FIXXED_G_OPEN     => :blue, 
    # FBA_Z_FIXXED_G_BOUNDED  => :blue,
)

MODEL_LABELS = Dict(
    :dyn => "\$ DYNAMIC \$",
    ME_FULL_POLYTOPE => "\$ ME^{\\langle 1,1 \\rangle} \$",
    ME_Z_EXPECTED_G_BOUNDED => "\$ ME^{\\langle 1,0 \\rangle} \$",
    ME_Z_FIXXED_G_BOUNDED => "\$ {ME^{\\langle 0,0 \\rangle}}_z \$", 
    ME_Z_OPEN_G_OPEN => "\$ ME^{\\langle 0,0 \\rangle} \$", 
    ME_Z_EXPECTED_G_EXPECTED => "\$ ME^{\\langle 2,0 \\rangle} \$", 
)

# :auto, :circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, 
# :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, 
# :star6, :star7, :star8, :vline, :hline, :+, :x
MODEL_MARKERS = Dict(
    :dyn => :circle,
    ME_FULL_POLYTOPE => :square,
    ME_Z_EXPECTED_G_BOUNDED => :utriangle,
    ME_Z_FIXXED_G_BOUNDED => :star, 
    ME_Z_OPEN_G_OPEN => :hex,
    ME_Z_EXPECTED_G_EXPECTED => :dtriangle,
)

MOD_LS = Dict(
    (ME_MODELS .=> :dash)...,
    (FBA_MODELS .=> :dot)...,
)

ES_COLORS = let
    ϵs = MINDEX[:ϵs]
    paltt = ColorSchemes.thermal
    colors = get.([paltt], ϵs / maximum(ϵs))
    Dict(ϵ => c for (ϵ, c) in zip(ϵs, colors))
end

# ----------------------------------------------------------------------------
# PLOTS
mysavefig(p, pname; params...) = 
    Dyn.mysavefig(p, pname, Dyn.plotsdir(), fileid; params...)

# ----------------------------------------------------------------------------
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

function plot_marginals!(p, MODsyms, rxn, Vl, D, ϵ, τ; 
        draw_dyn = true,
        draw_av = true, 
        draw_va = true, 
        sparams...
    )

    isdelta(M) = length(findall(!isequal(minimum(values(M))), M)) == 1

    # std fill
    if draw_dyn
        DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
        plot!(p, DyMs[rxn]; label = "", alpha = 0.0)
        if draw_va
            av = Dyn.av(DyMs[rxn])
            std = sqrt(Dyn.va(DyMs[rxn]))
            vspan!(p, [av - std, av + std], label = "",
                linecolor = :grey, fillcolor = :grey,
                alpha = 0.3
            )
        end
    end
    
    for MODsym in MODsyms
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        plot!(p, Ms[rxn]; label = "", alpha = 0.0)
        if draw_va && !isdelta(Ms[rxn])
            av = Dyn.av(Ms[rxn])
            std = sqrt(Dyn.va(Ms[rxn]))
            vspan!(p, [av - std, av + std], label = "",
                linecolor = color, fillcolor = color,
                alpha = 0.2
            )
        end
    end

    # means and std lines
    if draw_dyn
        av = Dyn.av(DyMs[rxn])
        std = sqrt(Dyn.va(DyMs[rxn]))
        draw_va && vline!(p, [av - std, av + std]; 
            label = "", sparams..., alpha = 0.5,
            color = :black, ls = :solid, lw = 3
        )
        draw_av && vline!(p, [av]; label = "", sparams..., 
            color = :black, ls = :solid, lw = 5
        )
    end
    
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        if !isdelta(Ms[rxn])
            av = Dyn.av(Ms[rxn])
            std = sqrt(Dyn.va(Ms[rxn]))
            draw_va && vline!(p, [av - std, av + std]; 
                label = "", color, sparams..., 
                ls, lw = 3, alpha = 0.3
            )
            draw_av && vline!(p, [av]; label = "", color, sparams..., 
                ls, lw = 5
            )
        end
    end

    # marginals
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        
        if isdelta(Ms[rxn])
            mx = findfirst(isone, Ms[rxn])
            xlb, xub = (keys(Ms[rxn]),) .|> (minimum, maximum)
            ylb, yub = (values(Ms[rxn]),) .|> (minimum, maximum)

            plot!(p, [xlb, xub], [ylb, ylb]; label = "", color, ls, sparams...)
            plot!(p, [mx, mx], [ylb, yub]; label = "", color, ls, sparams...)
        else
            plot!(p, Ms[rxn]; label = "", color, ls, sparams...)
        end
    end
    draw_dyn && plot!(p, DyMs[rxn]; 
        label = "", color = :black, sparams..., 
    )

    return p
end
plot_marginals!(p, MODsyms::Symbol, rxn, Vl, D, ϵ, τ; kwargs...) = 
    plot_marginals!(p, [MODsyms], rxn, Vl, D, ϵ, τ; kwargs...)
