module Chemostat_InSilico

    export PROJECT_NAME, PROJECT_DIR, FIGURES_DIR, DATA_DIR, CACHE_DIR
    export CH3_DATA_DIR, CH3_FIGURES_DIR

    include("Meta/Meta.jl")
    include("Utilities/Utilities.jl")
    # include("Polytopes/Polytopes.jl")
    # include("MonteCarlo/MonteCarlo.jl")
    # include("MaxEnt/MaxEnt.jl")
    include("LP_Implement/LP_Implement.jl")

    function __init__()
        _create_dirs()
    end

end
