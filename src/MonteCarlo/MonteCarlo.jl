module MonteCarlo

using Plots
using ProgressMeter
using Base.Threads
using ..Polytopes
using ..Utilities

export Cell, is_inpolytope, vg, vatp
export generate_random_cells, pick_cells, runMC

include("Cells.jl")
include("sampling_polytope.jl")
include("runMC.jl")

end