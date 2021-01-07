const PROJECT_NAME = "Chemostat_InSilico"
const PROJECT_DIR = abspath(joinpath(@__DIR__, "../.."))
const FIGURES_DIR = joinpath(PROJECT_DIR, "figures")
    const CH2_FIGURES_DIR = joinpath(FIGURES_DIR, "MCSim")
    const DYN_FIGURES_DIR = joinpath(FIGURES_DIR, "DynSim")
const DATA_DIR = joinpath(PROJECT_DIR, "data")
    const CH2_DATA_DIR = joinpath(DATA_DIR, "MCSim")
    const DYN_DATA_DIR = joinpath(DATA_DIR, "DynSim")
        const DYN_DATA_BUNDLE_FILE = joinpath(DYN_DATA_DIR, "dym_dat_bundle.bson")

const CACHE_DIR = joinpath(DATA_DIR, "cache")


function _create_dirs()
    for dir in [
            FIGURES_DIR, DATA_DIR, CACHE_DIR,
            CH2_DATA_DIR, CH2_FIGURES_DIR,
            DYN_DATA_DIR, DYN_FIGURES_DIR
        ]
        mkpath(dir)
    end
end