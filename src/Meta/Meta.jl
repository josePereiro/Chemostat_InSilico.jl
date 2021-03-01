# const PROJECT_NAME = "Chemostat_InSilico"
#     const DYN_PROJECT_NAME = "DynSim"
#     const MC_PROJECT_NAME = "MCSim"

# const PROJECT_DIR = abspath(joinpath(@__DIR__, "../.."))

#     const FIGURES_DIR = joinpath(PROJECT_DIR, "figures")
#         const CH2_FIGURES_DIR = joinpath(FIGURES_DIR, MC_PROJECT_NAME)
#         const DYN_FIGURES_DIR = joinpath(FIGURES_DIR, DYN_PROJECT_NAME)

#     const DATA_DIR = joinpath(PROJECT_DIR, "data")
#         const CH2_DATA_DIR = joinpath(DATA_DIR, MC_PROJECT_NAME)
#         const DYN_DATA_DIR = joinpath(DATA_DIR, DYN_PROJECT_NAME)
#             const DYN_DATA_INDEX_FILE = joinpath(DYN_DATA_DIR, "dym_dat_index.bson")
#             const MARGINALS_INDEX_FILE = joinpath(DYN_DATA_DIR, "marginals_dat_index.bson")

#     const CACHE_DIR = joinpath(DATA_DIR, "cache")
#         const DYN_CACHE_DIR = joinpath(CACHE_DIR, DYN_PROJECT_NAME)

# function _create_dirs()
#     for dir in [
#             FIGURES_DIR, DATA_DIR, CACHE_DIR,
#             CH2_DATA_DIR, CH2_FIGURES_DIR,
#             DYN_DATA_DIR, DYN_FIGURES_DIR, DYN_CACHE_DIR
#         ]
#         mkpath(dir)
#     end
# end