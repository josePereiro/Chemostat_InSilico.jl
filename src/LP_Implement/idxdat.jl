# TODO: package this
const IDXDAT_CACHE = Dict()
function idxdat(INDEX, dk, indexks...; cache = true)
    FILE = INDEX[:DFILE, indexks...]
    if FILE isa ITERABLE
        dat = []
        for F in FILE
            datum = deserialize(F)[dk...]
            push!(dat, datum)
        end
        return dat
    else
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