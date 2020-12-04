module Polytopes

using Plots
using Base.Threads
using ProgressMeter

export Polytope, ΔPolytope, AbstractPolytope
export setX!, getX, vg_range, vatp_range, vg_nbins, vatp_nbins, nbins, addXs!
export is_inpolytope, plot_polytope!, integrate_vatp_slicing
export vgL, vgU, vatpL, vatpU, vl, vr, z, vlU
export Δvg, Δvatp
export rand_vg, rand_vatp

include("AbstractPolytope.jl")
include("Polytope.jl")
include("ΔPolytope.jl")
include("plot_polytope.jl")
include("intagrate_vatp_splicing.jl")
include("rand.jl")

end