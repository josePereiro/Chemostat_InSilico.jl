module MonteCarlo

using Plots
using ProgressMeter
using Base.Threads
using ..Polytopes
import UtilsJL: get_chuncks

include("Cells.jl")
include("sampling_polytope.jl")
include("runMC.jl")

export Cell, is_inpolytope, vg, vatp
export generate_random_cells, pick_cells, runMC, generate_random_point, generate_similar_cell

end