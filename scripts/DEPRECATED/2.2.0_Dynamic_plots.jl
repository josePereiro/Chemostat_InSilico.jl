import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic

    using ProgressMeter
    using Plots
    using Plots.PlotMeasures
    using InteractiveUtils

    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Serialization
    using Statistics
    using Dates
    using Random
    using FileIO
    using ColorSchemes
end

# ----------------------------------------------------------------------------
# Meta
fileid = "2.2"
Base.first(v::Vector, i) = v[firstindex(v):(min(lastindex(v), i))]
Base.last(v::Vector, i) = v[max(firstindex(v), length(v) - i + 1):lastindex(v)]

# ----------------------------------------------------------------------------
# Load DAT
# INDEX[:DFILE, Vl, D, ϵ, τ]
INDEX = UJL.load_data(Dyn.procdir("dym_dat_index.bson"))
idxdat(dk, indexks...; cache = true) = Dyn.idxdat(INDEX, dk, indexks...; cache)
EXP_PARAMS = Iterators.product(INDEX[[:Vls, :Ds, :ϵs, :τs]]...);

# ----------------------------------------------------------------------------
# PLOTS
mysavefig(p, pname; params...) = 
    Dyn.mysavefig(p, pname, Dyn.plotsdir(), fileid; params...)

# ----------------------------------------------------------------------------
ϵs_colors = let
    ϵs = INDEX[:ϵs]
    paltt = ColorSchemes.thermal
    colors = get.([paltt], ϵs / maximum(ϵs))
    Dict(ϵ => c for (ϵ, c) in zip(ϵs, colors))
end

