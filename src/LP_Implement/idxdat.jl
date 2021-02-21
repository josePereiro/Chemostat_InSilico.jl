# TODO: package this
const IDXDAT_CACHE = Dict()
function idxdat(INDEX, dk::Vector, indexks...; cache = true, emptycache = false)
    emptycache && empty!(IDXDAT_CACHE)
    FILE = INDEX[:DFILE, indexks...]
    if FILE isa ITERABLE
        dat = []
        for F in FILE
            F = joinpath(InCh.PROJECT_DIR, F)
            datum = deserialize(F)[dk...]
            push!(dat, datum)
        end
        return dat
    else
        FILE = joinpath(InCh.PROJECT_DIR, FILE)
        iscached = cache && FILE == get(IDXDAT_CACHE, :FCACHED, nothing)
        return iscached ?
            IDXDAT_CACHE[:DATCACHED][dk...] :
        begin
            dat = deserialize(FILE)
            IDXDAT_CACHE[:FCACHED] = FILE
            IDXDAT_CACHE[:DATCACHED] = dat
            dat[dk...]
        end
    end
end