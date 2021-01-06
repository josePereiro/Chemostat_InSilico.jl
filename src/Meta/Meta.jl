const PROJECT_NAME = "Chemostat_InSilico"
const PROJECT_DIR = abspath(joinpath(@__DIR__, "../.."))
const FIGURES_DIR = joinpath(PROJECT_DIR, "figures")
    const CH2_FIGURES_DIR = joinpath(FIGURES_DIR, "Ch2Sim")
    const CH3_FIGURES_DIR = joinpath(FIGURES_DIR, "Ch3Sim")
const DATA_DIR = joinpath(PROJECT_DIR, "data")
    const CH2_DATA_DIR = joinpath(DATA_DIR, "Ch2Sim")
    const CH3_DATA_DIR = joinpath(DATA_DIR, "Ch3Sim")
        const CH3_DAT_BUNDLE_FILE = joinpath(CH3_DATA_DIR, "dat_bundle.bson")

const CACHE_DIR = joinpath(DATA_DIR, "cache")


function _create_dirs()
    for dir in [
            FIGURES_DIR, DATA_DIR, CACHE_DIR,
            CH2_DATA_DIR, CH2_FIGURES_DIR,
            CH3_DATA_DIR, CH3_FIGURES_DIR
        ]
        mkpath(dir)
    end
end