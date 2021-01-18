# TODO: package this
FCACHED = nothing
DATCACHED = nothing
function idxdat(INDEX, dk, indexks...)
    FILE = INDEX[:DFILE, indexks...]
    if FILE isa ITERABLE
        dat = []
        for F in FILE
            datum = deserialize(F)[dk]
            push!(dat, datum)
        end
        return dat
    else
        if FILE == FCACHED
            return DATCACHED[dk]
        else
            global FCACHED = FILE
            dat = deserialize(FILE)
            global DATCACHED = dat
            return dat[dk]
        end
    end
end