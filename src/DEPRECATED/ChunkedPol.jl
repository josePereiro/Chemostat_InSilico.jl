struct ChunkedPol{T}
    net::MetNet
    box::BoxGrid
    chnklen::Int
    len::Int
    chnk_files::Vector{String}

    function ChunkedPol(net, frees, δ::Int; 
            chnkthMB = 500 , # chunk threshold in bytes
            tolf = 0.01, 
            dir = procdir(), 
            nthrs = nthreads(),
            elT = Float16,
            ignorecache = false
        )

        # chnklen
        chnkthB = chnkthMB * (1024 * 1024)
        Bpu = length(frees) * sizeof(elT)
        chnklen = div(chnkthB, Bpu)

        # check cache
        pol_cfile = _pol_chunk_sfile(net, chnkthMB, tolf; dir)
        (isfile(pol_cfile) && !ignorecache) && return ldat(pol_cfile; verbose = false)

        freeis = rxnindex(net, frees)
        box = BoxGrid(net, freeis, δ)
        tol = Δv(net, freeis) .* tolf
        chnk_files = String[]

        boxlen = length(box)
        (chnklen > boxlen) && (chnklen = boxlen)

        @info("Setup", chnklen, boxlen, boxlen/chnklen, chnkthMB, nthrs)

        chnk_no = 0
        wl = ReentrantLock()
        ch = UtilsJL.GeneralUtils.chunkedChannel(box; 
            chnklen = max(10_000, div(chnklen, nthrs)), 
            buffsize = nthrs
        )
        
        fea_count = 0
        fea_count_estimate = 0
        
        prog = Progress(boxlen; dt = 0.5, desc = "Collecting Pol... ")
        showvalues() = [
            ("chunks", length(chnk_files)), 
            ("fea_count_estimate", fea_count_estimate), 
            ("boxlen", boxlen), 
            ("polV/boxV", fea_count_estimate/boxlen),
            ("thid", threadid()),
        ]

        @threads for _ in 1:(nthrs + 1)
            thid = threadid()
            (thid == 1) && continue

            pol_chnk = Vector{Vector{elT}}()
            idx = 0
            net_th = deepcopy(net)
            
            for box_chnk in ch
                isempty(pol_chnk) && 
                    (pol_chnk = Vector{Vector{elT}}(undef, chnklen))
                
                for f64v in box_chnk
                    if is_feasible(net_th, freeis, f64v; tol)
                        idx += 1
                        pol_chnk[idx] = collect(elT, f64v)
                        
                        # end of pol_chnk
                        if (idx == chnklen)
                            lock(wl) do
                                chnk_file = _save_pol_chunk(pol_chnk, pol_cfile, chnk_no; dir)
                                push!(chnk_files, chnk_file)
                                fea_count += chnklen
                                fea_count_estimate = fea_count
                                @info("Done!", thid, chnklen, idx, chnk_no, fea_count)
                                chnk_no += 1
                            end
                            idx = 0
                        end
                        
                        fea_count_estimate += 1
                    end # if is_feasible

                    next!(prog; showvalues)

                end # for f64v in box_chnk
            end # for box_chnk in ch

            # last chunk
            if idx != 0
                resize!(pol_chnk, idx)
                lock(wl) do
                    chnk_file = _save_pol_chunk(pol_chnk, pol_cfile, chnk_no)
                    push!(chnk_files, chnk_file)
                    fea_count += idx
                    fea_count_estimate = fea_count
                    # @info("Done!", thid, chnklen, idx, chnk_no, fea_count)
                    chnk_no += 1
                end
            end

        end # for _ in _threads

        finish!(prog)

        chnkpol = new{elT}(
            net, box, chnklen, fea_count, 
            chnk_files
        )
        sdat(chnkpol, pol_cfile; verbose = false)
        return chnkpol
    end
end


_pol_chunk_sfile(net, chnklen, tolf; dir) = 
    Ass.dfname([dir], "pol_chunk", (;net = hash(net), chnklen, tolf), ".jls")

_save_pol_chunk(pol_chnk, pol_cfile, chnk_no; dir = procdir()) =
    Ass.sdat(pol_chnk, [dir], "pol_chunk", chnk_no, (;hash = hash(pol_cfile)), ".jls"; verbose = false)


## ------------------------------------------------------
function Base.show(io::IO, ckp::ChunkedPol{T}) where {T}
    println(io, "ChunkedPol{", T, "}")
    println(io, "frees: [", join(string.(ckp.box.dims_ids), ", "), "]")
    println(io, "chunks: ", length(ckp.chnk_files))
    println(io, "chunklen: ", ckp.chnklen)
    println(io, "len: ", ckp.len)
    println(io, ckp.net)
end
Base.show(ckp::ChunkedPol) = show(stdout, ckp)

## ------------------------------------------------------
# Iterator Interface
Base.IteratorSize(::ChunkedPol) = Base.HasLength()
Base.IteratorEltype(::ChunkedPol{T}) where {T} = Base.HasEltype()
Base.eltype(ckp::ChunkedPol{T}) where {T} = Vector{T}
Base.length(ckp::ChunkedPol) = ckp.len
Base.size(ckp::ChunkedPol) = tuple(ckp.len)
Base.size(ckp::ChunkedPol, dim) = size(ckp)[dim]

function _load_pol_chunk(ckp::ChunkedPol{T}, file)::Vector{Vector{T}} where {T}
    !isfile(file) && error("Cache file missing, ", file)
    dat = Ass.ldat(file; verbose = false)
    isempty(dat) && error("Cache is empty")
    return dat
end

## ------------------------------------------------------
function Base.iterate(ckp::ChunkedPol)
    
    # dat file
    chnk_files = ckp.chnk_files
    isempty(chnk_files) && error("chnk_files is empty")
    file = first(chnk_files)
    next_file_idx = 2
    
    # load dat
    dat = _load_pol_chunk(ckp, file)
    next_dat_idx = 2
    datlen = length(dat)
    elem = first(dat)

    return (elem, (dat, datlen, next_dat_idx, file, next_file_idx))
end

function Base.iterate(ckp::ChunkedPol, state)
    (dat, datlen, next_dat_idx, file, next_file_idx) = state
    
    # next elem
    if next_dat_idx <= datlen
        elem = dat[next_dat_idx]
        next_dat_idx += 1
        return (elem, (dat, datlen, next_dat_idx, file, next_file_idx))
    end
    
    # load new dat
    chnk_files = ckp.chnk_files
    if next_file_idx <= length(chnk_files)
        
        # dat file
        file = chnk_files[next_file_idx]
        next_file_idx += 1
        
        # load dat
        dat = _load_pol_chunk(ckp, file)
        datlen = length(dat)
        elem = first(dat)
        next_dat_idx = 2

        return (elem, (dat, datlen, next_dat_idx, file, next_file_idx))
    end

    # end iter
    return nothing

end
