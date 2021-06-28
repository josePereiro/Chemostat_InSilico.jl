function poly_vol_board(M, Ds, cgD_Xs; 
        bins = 10, nths = nthreads(), 
        freeth::Int = 0, 
        clear_cache = true
    )

    # Check cache
    chash = (hash(M.net), hash(Ds), hash(cgD_Xs), bins)
    clear_cache && UJL.delete_cache(chash)
    cdat = UJL.load_cache(chash; verbose = false)
    !isnothing(cdat) && return cdat

    psize = zeros(length(Ds), length(cgD_Xs))

    # feeding task
    Ch = Channel(nths) do Ch_
        @showprogress for (Di, D) in enumerate(Ds)
            for (cgD_Xi, cgD_X) in enumerate(cgD_Xs)
                put!(Ch_, (Di, D, cgD_Xi, cgD_X))
            end
        end 
    end

    @threads for _ in 1:nths
        thid = threadid()
        nths > 1 && thid == freeth && continue # Let it free

        net = deepcopy(M.net)
        for (Di, D, cgD_Xi, cgD_X) in Ch
            psize[Di, cgD_Xi] = poly_vol(M, D, cgD_X; net)
        end
    end

    # cache
    UJL.save_cache(chash, psize)
    
    return psize
end

function poly_vol(M, D, cgD_X; 
        net = M.net
    )

    biom_idx = M.obj_idx
    vg_idx = M.vg_idx

    vgL, vgU = fixxing(net, biom_idx, D) do
        fva(net, vg_idx)
    end
    ifeasible = cgD_X >= vgL && (vgL != vgU && vgL != 0.0)
    ifeasible ? cgD_X - vgL : NaN
end