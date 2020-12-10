module Utilities

using ProgressMeter
using DrWatson
using Printf

export get_chuncks, grad_desc
export save_data, load_data, mysavename

include("get_chuncks.jl")
include("gradient_descent.jl")
include("save_and_load.jl")

end

