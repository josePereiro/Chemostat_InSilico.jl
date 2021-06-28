function proj2d(M, ider1, ider2; bins = 100)
    net = deepcopy(M.net)
    idx1 = rxnindex(net, ider1)
    idx2 = rxnindex(net, ider2)

    proj2d = Dict()
    idx1L, idx1U = fva(net, idx1)
    range1 = range(idx1L, idx1U; length = bins)
    for val1 in range1
        net.lb[idx1] = net.ub[idx1] = val1
        idx2L, idx2U = fva(net, idx2)
        proj2d[val1] = (idx2L, idx2U)
    end
    proj2d
end