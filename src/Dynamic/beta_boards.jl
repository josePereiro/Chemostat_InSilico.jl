## -----------------------------------------------------------------------------------------------
#=
Produce z_beta vs _vg_beta boards that can be use to
plot head maps of entropy(H) or vg and z averages
=#
function get_beta_boards(M, zbetas, vgbetas; 
        nths = nthreads(), 
        freeth = 0,
        clear_cache = true
    )

    # Check cache
    chash = (hash(M.net), hash(zbetas), hash(vgbetas))
    clear_cache && UJL.delete_cache(chash)
    cdat = UJL.load_cache(chash;verbose = false)
    !isnothing(cdat) && return cdat

    obj_idx = M.obj_idx
    LP_cache = vgvatp_cache(M)
    z(vatp, vg) = LP_cache[vatp][vg][obj_idx]
    vg(vatp, vg) = vg

    zN, vgN = length(zbetas), length(vgbetas)

    H_board = zeros(zN, vgN)
    vg_board = similar(H_board)
    z_board = similar(H_board)
    prog = Progress(length(H_board))

    # feed task
    Ch = Channel(nthreads()) do Ch_
        for (z_bi, z_beta) in enumerate(zbetas)
            for (vg_bi, vg_beta) in enumerate(vgbetas)
                put!(Ch_, (z_bi, z_beta, vg_bi, vg_beta))
            end
        end
    end

    # Compute boards
    @threads for _ in 1:nths
        thid = threadid()
        nths > 1 && thid == freeth && continue # Let it free
        
        M0 = deepcopy(M)
        PME = get_join(M0)
        for (z_bi, z_beta, vg_bi, vg_beta) in Ch
            
            try;
                PME = get_join!(M0, PME) do _vatp, _vg
                    exp(z_beta * z(_vatp, _vg) + vg_beta * vg(_vatp, _vg))
                end
                H = entropy(PME)

                H_board[z_bi, vg_bi] = H
                vg_board[z_bi, vg_bi] = ave_over(PME) do vatp, vg
                    vg
                end
                z_board[z_bi, vg_bi] = ave_over(PME) do vatp, vg
                    z(vatp, vg)
                end

                next!(prog; showvalues = [
                        (:z_beta, z_beta),    
                        (:vg_beta, vg_beta),    
                        (:H, H),    
                        (:thid, thid),    
                    ]
                )
            catch err
                @warn("Error at", z_beta, vg_beta)
                rethrow(err)
            end
        end

    end # for thid

    # cache
    UJL.save_cache(chash, (H_board, z_board, vg_board))

    return H_board, z_board, vg_board
end

getbeta(order) = sign(order) * 10.0 ^ abs(order)