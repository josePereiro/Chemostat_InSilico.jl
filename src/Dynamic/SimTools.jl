## ------------------------------------------------------
# globals
const UNDONE_SIM_STATUS = :UNDONE
const DEAD_SIM_STATUS = :DEAD
const EXPLODED_SIM_STATUS = :EXPLODED
const NITERS_SIM_STATUS = :NITERS
const STST_SIM_STATUS = :STST
const STATUS_FILE_EXT = ".status.jls"
const BATCH_FILE_EXT = ".batch.jls"
const SIM_FILE_EXT = ".simdat.jls"
const BFILES_FILE_EXT = ".bfiles.jls"
const ME_FILE_EXT = ".maxent.jls"

## ------------------------------------------------------
# utils
status_file(simid, simparams) = procdir(simid, simparams, STATUS_FILE_EXT)
set_status(status, simid, simparams) = (sdat(status, status_file(simid, simparams); verbose = false); status)
get_status(simid, simparams) = ldat(() -> UNDONE_SIM_STATUS, status_file(simid, simparams); verbose = false)

## ------------------------------------------------------
# batch files
bfiles_file(simid, simparams) = procdir(simid, simparams, BFILES_FILE_EXT)
save_bfiles(bfiles, simid, simparams) = sdat(bfiles, bfiles_file(simid, simparams); verbose = false)
load_bfiles(simid, simparams) = ldat(bfiles_file(simid, simparams); verbose = false)

## ------------------------------------------------------
# batch
batch_file(simid, simparams, it) = procdir(simid, simparams, (;it), BATCH_FILE_EXT)
save_batch(dat, simid, simparams, it) = 
   (bfile = batch_file(simid, simparams, it); sdat(dat, bfile; verbose = false); bfile)
load_batch(simid, simparams, it) = ldat(batch_file(simid, simparams, it); verbose = false)

## ------------------------------------------------------
# simdat
simdat_file(simid, simparams) = procdir(simid, simparams, SIM_FILE_EXT)
save_simdat(S, simid, simparams) = sdat(S, simdat_file(simid, simparams); verbose = false)
load_simdat(simid, simparams) = ldat(simdat_file(simid, simparams); verbose = false)

## ------------------------------------------------------
# maxent
maxent_file(simid, MEmode, simparams) = 
   procdir(simid, simparams, (;MEmode = string(MEmode)), ME_FILE_EXT)
save_maxent(dat, simid, MEmode, simparams) = 
   sdat(dat, maxent_file(simid, MEmode, simparams); verbose = false)
load_maxent(simid, MEmode, simparams) = 
   ldat(maxent_file(simid, MEmode, simparams); verbose = false)

## ------------------------------------------------------
function inparams(nparams::Dict; qparams...)
   qparams = Dict(qparams...)
   isempty(qparams) && return false
   isempty(nparams) && return false
   for (k, qparam) in qparams
      (get(nparams, string(k), nothing) != qparam) && return false
   end
   return true
end
 
## ------------------------------------------------------
function collect_ts(tsid, simid, qparams)
   dat_pool = Dict()
   bfiles = load_bfiles(simid, qparams)
   t0 = 1
   for bname in bfiles
      # batch
      batch = ldat(bname; verbose = false)
      bs = getfield(batch, tsid)
      @assert (bs isa Vector)
      blen = length(bs)
      
      bt = (t0:(t0 + blen - 1)) .* batch.push_frec
      ts = get!(dat_pool, :ts, Float64[])
      push!(ts, bt...)
      bdat = identity.(bs)
      ds = get!(dat_pool, :ds, eltype(bdat)[])
      push!(ds, bdat...)
      t0 += length(bs) + 1
   end
   ts, ds = get(dat_pool, :ts, []), get(dat_pool, :ds, [])
   return (ts, ds)
end