## ------------------------------------------------------
function move_frees!(net, frees::Vector)
    idxs = Dyn.rxnindex.([net], frees)
    M, N = size(net)
    offset = N - rank(net.S)
    (length(frees) != offset) && 
        error(
            "The system has ", offset, 
            " degree of freedom, provided only ", 
            length(frees), " frees"
        )
    
    for (i, oldidx) in enumerate(idxs)
        newidx = (N - offset) + i

        # S
        Sbk = net.S[:, newidx]
        net.S[:, newidx] .= net.S[:, oldidx]
        net.S[:, oldidx] .= Sbk

        # vectors
        for field in [:lb, :ub, :c, :rxns]
            v = getfield(net, field)
            bk = v[newidx]
            v[newidx] = v[oldidx]
            v[newidx] = v[oldidx]
            v[oldidx] = bk
        end
    end

    return net
end

## ------------------------------------------------------
function reduced_net(net, frees)
    net = deepcopy(net)
    move_frees!(net, frees)
    frees = Dyn.rxnindex.([net], frees)

    # echelon
    F = lu([net.S net.b])
    echSb = rref(F.L * F.U)
    return Dyn.MetNet(;
        S = echSb[:, frees], 
        mets = net.mets[F.p],
        b = echSb[:, end],
        rxns = net.rxns[frees],
        lb = net.lb[frees],
        ub = net.ub[frees],
        c = net.c[frees], 
    )
end