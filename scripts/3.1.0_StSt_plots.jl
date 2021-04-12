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

# ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
minmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

# ----------------------------------------------------------------------------
# Prepare network

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

ALL_MODELS = [
    # ME_Z_OPEN_G_OPEN, ME_Z_OPEN_G_BOUNDED, 
    # ME_Z_EXPECTED_G_OPEN, 
    # ME_Z_EXPECTED_G_BOUNDED, 
    ME_FULL_POLYTOPE,
    # ME_Z_EXPECTED_G_MOVING,
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
    ME_FULL_POLYTOPE         => :brown,
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
    ME_FULL_POLYTOPE         => :dash,
    ME_Z_EXPECTED_G_MOVING => :dash,
    ME_Z_FIXXED_G_OPEN      => :dash, 
    ME_Z_FIXXED_G_BOUNDED   => :dash, 
    ME_Z_EXPECTED_G_EXPECTED   => :dash, 

    FBA_Z_OPEN_G_OPEN       => :dot, 
    FBA_Z_OPEN_G_BOUNDED    => :dot, 
    FBA_Z_FIXXED_G_OPEN     => :dot, 
    FBA_Z_FIXXED_G_BOUNDED  => :dot,
)
 
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

function plot_marginals!(p, MODsyms, rxn, Vl, D, ϵ, τ; draw_av = true, sparams...)

    # std fill
    DyMs = idxdat([:DyMs], Vl, D, ϵ, τ)
    plot!(p, DyMs[rxn]; label = "", alpha = 0.0)
    if draw_av
        av = Dyn.av(DyMs[rxn])
        std = sqrt(Dyn.va(DyMs[rxn]))
        vspan!(p, [av - std, av + std], label = "",
            linecolor = :grey, fillcolor = :grey,
            alpha = 0.3
        )
    end
    
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        plot!(p, Ms[rxn]; label = "", alpha = 0.0)
        if draw_av
            av = Dyn.av(Ms[rxn])
            std = sqrt(Dyn.va(Ms[rxn]))
            vspan!(p, [av - std, av + std], label = "",
                linecolor = color, fillcolor = color,
                alpha = 0.3
            )
        end
    end

    # means and std lines
    if draw_av
        av = Dyn.av(DyMs[rxn])
        std = sqrt(Dyn.va(DyMs[rxn]))
        vline!(p, [av - std, av + std]; 
            label = "", sparams..., 
            color = :black, ls = :solid, lw = 3
        )
        vline!(p, [av]; label = "", sparams..., 
            color = :black, ls = :solid, lw = 5
        )
    end
    
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        if draw_av
            av = Dyn.av(Ms[rxn])
            std = sqrt(Dyn.va(Ms[rxn]))
            vline!(p, [av - std, av + std]; 
                label = "", color, sparams..., 
                ls, lw = 3
            )
            vline!(p, [av]; label = "", color, sparams..., 
                ls, lw = 5
            )
        end
    end

    # marginals
    for MODsym in MODsyms
        ls = MOD_LS[MODsym] 
        color = MOD_COLORS[MODsym]
        Ms = idxdat([MODsym, :Ms], Vl, D, ϵ, τ)
        plot!(p, Ms[rxn]; label = "", color, ls, sparams...)
    end
    plot!(p, DyMs[rxn]; label = "", sparams..., color = :black)

    return p
end
plot_marginals!(p, MODsyms::Symbol, rxn, Vl, D, ϵ, τ; sparams...) = 
    plot_marginals!(p, [MODsyms], rxn, Vl, D, ϵ, τ; sparams...)


