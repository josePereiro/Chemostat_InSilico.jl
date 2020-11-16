module Polytopes

using Plots
using Base.Threads
using ProgressMeter

export Polytope, is_inpolytope, plot_polytope, integrate_vatp_slicing
export vgL, vgU, vatpL, vatpU
export vgL, vgU, vatpL, vatpU
export Δvg, Δvatp
export rand_vg, rand_vatp

include("type.jl")
include("polytope_bounds.jl")
include("plot_polytope.jl")
include("intagrate_vatp_splicing.jl")
include("rand.jl")

end