module Chemostat_InSilico

    import UtilsJL
    UtilsJL.gen_top_proj(@__MODULE__)

    # include("Polytopes/Polytopes.jl")
    # include("MaxEnt/MaxEnt.jl")
    # include("MonteCarlo/MonteCarlo.jl")
    include("Dynamic/Dynamic.jl")

    function __init__()
        UtilsJL.create_proj_dirs(@__MODULE__)
    end

end
