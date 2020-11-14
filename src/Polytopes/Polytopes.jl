module Polytopes

using Plots
using Base.Threads
using ProgressMeter

export Polytope, is_inpolytope, plot_polytope, integrate_vatp_slicing
export vg_local_min, vg_local_max, vatp_local_min, vatp_local_max
export vg_global_min, vg_global_max, vatp_global_min, vatp_global_max
export Δvg, Δvatp

include("type.jl")
include("polytope_bounds.jl")
include("plot_polytope.jl")
include("intagrate_vatp_splicing.jl")

end