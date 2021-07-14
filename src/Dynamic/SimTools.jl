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
 function filter_dfnames(f::Function, names::Vector)
    
    return filter(names) do name
       heads, nparams, ext = Ass.parse_dfname(name)
       f(heads, nparams, ext)
    end
 end
 
 ## ------------------------------------------------------
function _to_symbol_dict(d::Dict{KT, VT}) where {KT, VT}
    sd = Dict{Symbol, VT}()
    for (kstr, dat) in d
        sd[Symbol(kstr)] = dat
    end
    sd
end

## ------------------------------------------------------
 function read_batch_files(simid, qparams)
    names = readdir(procdir())

    dfname_ = Ass.dfname("_test", qparams)
    qparams_ = _to_symbol_dict(Ass.dfparams(dfname_))

    names = filter_dfnames(names) do heads, params, ext
       (ext != BATCH_FILE_EXT) && return false
       (first(heads) != simid) && return false
       !inparams(params; qparams_...) && return false
       return true
    end
    sort!(names; by = (n) -> Ass.dfparam(n, "it", -Inf))
 end
 
 ## ------------------------------------------------------
 function collect_ts(tsid, simid, qparams)
    dat_pool = Dict()
    bfiles = read_batch_files(simid, qparams)
    t0 = 1
    for bname in bfiles
       # batch
       batch = ldat(bname; verbose = false)
       bs = getfield(batch, tsid)
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