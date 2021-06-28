module Dynamic

    using Core: Vector
using Plots
    import UtilsJL
    const Ass = UtilsJL.ProjAssistant
    Ass.gen_sub_proj(@__MODULE__)
    import MathProgBase.HighLevelInterface: linprog
    import Clp: ClpSolver
    using ProgressMeter

    include("discretize.jl")
    include("BoxGrid.jl")
    include("LP.jl")
    include("MetNets.jl")
    include("ToyModelD3.jl")
    include("PlotUtils.jl")
    # include("echelon.jl") # WIP

    function __init__()
        Ass.create_proj_dirs(@__MODULE__)
    end

end