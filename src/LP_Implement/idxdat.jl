# TODO: package this
const IDXDAT_CACHE = Dict()
function idxdat(INDEX, dk::Vector, indexks...; cache = true, emptycache = false)
    
    # CACHE
    thid = threadid()
    emptycache && empty!(IDXDAT_CACHE)
    TCACHE = get!(IDXDAT_CACHE, thid, Dict())

    DFILE = INDEX[:DFILE, indexks...]
    if DFILE isa ITERABLE
        dat = []
        for F in DFILE
            F = joinpath(InCh.PROJECT_DIR, F)
            datum = try deserialize(F)
                catch err; @warn("Error", F); rethrow(err)
            end
            push!(dat, datum[dk...])
        end
        return dat
    else
        F = joinpath(InCh.PROJECT_DIR, DFILE)
        iscached = F == get(TCACHE, :FCACHED, nothing)
        if cache && iscached 
            return TCACHE[:DATCACHED][dk...]
        else
            dat = try deserialize(F)
                catch err; @warn("Error", F); rethrow(err)
            end
            TCACHE[:FCACHED] = F
            TCACHE[:DATCACHED] = dat
            return dat[dk...]
        end
    end
end

function idxdat(;emptycache = true, emptytcache = true)
    emptycache && empty!(IDXDAT_CACHE)
    emptytcache && empty!(get!(IDXDAT_CACHE, threadid(), Dict()))
end