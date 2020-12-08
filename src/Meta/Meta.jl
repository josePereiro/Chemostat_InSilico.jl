const PROJECT_NAME = "Chemostat_Dynamics"
const PROJECT_DIR = abspath(joinpath(@__DIR__, "../.."))
const FIGURES_DIR = joinpath(PROJECT_DIR, "figures")
const DATA_DIR = joinpath(PROJECT_DIR, "data")
const CACHE_DIR = joinpath(DATA_DIR, "cache")

    const CH3_DATA_DIR = joinpath(DATA_DIR, "Ch3Sim")
    const CH3_FIGURES_DIR = joinpath(FIGURES_DIR, "Ch3Sim")

function _create_dirs()
    for dir in [
            FIGURES_DIR, DATA_DIR, CACHE_DIR,
            CH3_DATA_DIR, CH3_FIGURES_DIR
        ]
        mkpath(dir)
    end
end