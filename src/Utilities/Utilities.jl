module Utilities

using ProgressMeter
export get_chuncks, grad_desc

include("get_chuncks.jl")
include("gradient_descent.jl")

end

