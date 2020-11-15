module Utilities

using ProgressMeter
using DrWatson

export get_chuncks, grad_desc
export save_data, load_data

include("get_chuncks.jl")
include("gradient_descent.jl")
include("save_and_load.jl")

end

