## ------------------------------------------------------
# globals
const UNDONE_SIM_STATUS = :UNDONE
const DEAD_SIM_STATUS = :DEAD
const EXPLODED_SIM_STATUS = :EXPLODED
const NITERS_SIM_STATUS = :NITERS
const STST_SIM_STATUS = :STST
const STATUS_FILE_EXT = ".status.jls"
const BATCH_FILE_EXT = ".batch.jls"
const SIM_FILE_EXT = ".sim.jls"

## ------------------------------------------------------
# utils
status_file(simid, simparams) = procdir(simid, simparams, STATUS_FILE_EXT)
set_status(status, simid, simparams) = (sdat(status, status_file(simid, simparams); verbose = false); status)
get_status(simid, simparams) = ldat(() -> UNDONE_SIM_STATUS, status_file(simid, simparams); verbose = false)

## ------------------------------------------------------
# save batch
batch_file(simid, simparams, it) = procdir(simid, simparams, (;it), BATCH_FILE_EXT)
save_batch(dat, simid, simparams, it) = sdat(dat, batch_file(simid, simparams, it); verbose = false)
load_batch(simid, simparams, it) = ldat(batch_file(simid, simparams, it); verbose = false)

## ------------------------------------------------------
# save sim
sim_file(simid, simparams) = procdir(simid, simparams, SIM_FILE_EXT)
save_sim(S, simid, simparams) = sdat(S, sim_file(simid, simparams); verbose = false)
load_sim(simid, simparams) = ldat(sim_file(simid, simparams); verbose = false)
