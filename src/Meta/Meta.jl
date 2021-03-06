const PROJECT_NAME = "Chemostat_Dynamics"
const PROJECT_DIR = abspath(joinpath(@__DIR__, "../.."))
const FIGURES_DIR = joinpath(PROJECT_DIR, "figures")
const DATA_DIR = joinpath(PROJECT_DIR, "data")

function _create_dirs()
    for dir in [FIGURES_DIR, DATA_DIR]
        mkpath(dir)
    end
end