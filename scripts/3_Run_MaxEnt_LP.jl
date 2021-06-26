import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

# ----------------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--new-dat"
        help = "ignore disk stored DAT"   
        action = :store_true
    "--redo-me"
        help = "recompute MaxEnt part"   
        action = :store_true
    "--redo-fba"
        help = "recompute LP part"   
        action = :store_true
    "--dry-run"
        help = "do all but change nothing on disk"   
        action = :store_true
end

if isinteractive()
    # Dev values
    REDO_MAXENT = false
    REDO_FBA = false
    DRY_RUN = false
else
    parsed_args = parse_args(set)
    REDO_MAXENT = parsed_args["redo-me"]
    REDO_FBA = parsed_args["redo-fba"]
    DRY_RUN = parsed_args["dry-run"]
end

# ----------------------------------------------------------------------------
@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic

    import UtilsJL
    const UJL = UtilsJL
    using Base.Threads
    using Dates
    using Serialization
    using Random
end

# Compat
# @eval InCh const LP_Implement = Dynamic

## ----------------------------------------------------------------------------
# Load and clear DAT
# DINDEX [Vl, D, ϵ, τ] 
DINDEX_FILE = Dyn.procdir("dym_dat_index.bson")
DINDEX = UJL.load_data(DINDEX_FILE) # Dynamic index
DATA_FILE_PREFFIX = "marginal_dat"
dat_file(;sim_params...) = 
    Dyn.procdir(UJL.mysavename(DATA_FILE_PREFFIX, "jls"; sim_params...)) 
dyn_dat(dk, indexks...; cache = false, emptycache = true) = 
    Dyn.idxdat(DINDEX, dk, indexks...; cache, emptycache)

# ----------------------------------------------------------------------------
const WLOCK = ReentrantLock()

# ----------------------------------------------------------------------------
# METHOD VARIANTS
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_OPEN_G_BOUNDED       = :ME_Z_OPEN_G_BOUNDED        # 

const ME_Z_EXPECTED_G_OPEN      = :ME_Z_EXPECTED_G_OPEN       # Match ME and Dy biom average
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED   # 
const ME_FULL_POLYTOPE          = :ME_FULL_POLYTOPE           # 
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug

const ME_Z_FIXXED_G_OPEN        = :ME_Z_FIXXED_G_OPEN         # Fix biom around observed
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

const FBA_Z_OPEN_G_OPEN       = :FBA_Z_OPEN_G_OPEN
const FBA_Z_OPEN_G_BOUNDED    = :FBA_Z_OPEN_G_BOUNDED
const FBA_Z_FIXXED_G_OPEN     = :FBA_Z_FIXXED_G_OPEN 
const FBA_Z_FIXXED_G_BOUNDED  = :FBA_Z_FIXXED_G_BOUNDED
const FBA_Z_FIXXED_MAX_G      = :FBA_Z_FIXXED_MAX_G
const FBA_Z_FIXXED_MIN_G      = :FBA_Z_FIXXED_MIN_G

## ----------------------------------------------------------------------------
# set up functions
include("3.0.1_run_ME.jl")
include("3.0.2_run_FBA.jl")

## ----------------------------------------------------------------------------
# run
include("3.0.3_compute_and_save.jl")

## ----------------------------------------------------------------------------
# SAVING
MINDEX_FILE = Dyn.procdir("marg_dat_index.bson")
DRY_RUN || UJL.save_data(MINDEX_FILE, MINDEX)


