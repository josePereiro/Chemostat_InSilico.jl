module Chemostat_InSilico

    export PROJECT_NAME, PROJECT_DIR, FIGURES_DIR, DATA_DIR, CACHE_DIR
    export DYN_DATA_DIR, DYN_FIGURES_DIR

    include("Meta/Meta.jl")
    include("Polytopes/Polytopes.jl")
    include("MonteCarlo/MonteCarlo.jl")
    include("MaxEnt/MaxEnt.jl")
    include("LP_Implement/LP_Implement.jl")

    function __init__()
        _create_dirs()
    end

end
