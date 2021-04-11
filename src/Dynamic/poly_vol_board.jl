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

    biom_idx = M.obj_idx
    vg_idx = M.vg_idx

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
            exp_biom = D
            vgL, vgU = fixxing(net, biom_idx, exp_biom) do
                fva(net, vg_idx)
            end
            ifeasible = cgD_X >= vgL && (vgL != vgU && vgL != 0.0)
            psize[Di, cgD_Xi] = ifeasible ? cgD_X - vgL : NaN

        end
    end

    # cache
    UJL.save_cache(chash, psize)
    
    return psize
end