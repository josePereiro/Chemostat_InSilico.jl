function pos_defined(S, b, lb, ub, c, rxns, mets)
    M, N = size(S)
    bkwdis = findall(lb .< 0)
    bkwd_count = length(bkwdis)

    newN = N + bkwd_count
    newS = similar(S, M, newN) 
    newS[1:M, 1:N] .= S
    newb = copy(b)
    newmets = copy(mets)
    newlb = similar(lb, newN)
    newlb[1:N] .= lb
    newub = similar(ub, newN)
    newub[1:N] .= ub
    newrxns = similar(rxns, newN)
    newrxns[1:N] .= rxns
    newc = zeros(newN)
    newc[1:N] .= c
    # create bkwrs
    c = 1
    for fwdi in bkwdis
        bkwdi = N + c

        # bounds
        newub[bkwdi] = -lb[fwdi]
        newlb[fwdi] = newlb[bkwdi] = 0.0
        newS[:, bkwdi] .= -S[:, fwdi]
        newrxns[bkwdi] = string(rxns[fwdi], "_bkwd")
        c += 1
    end
    MetNet(newS, newb, newlb, newub, newc, newrxns, newmets)
end
pos_defined(net::MetNet) = pos_defined(net.S, net.b, net.lb, net.ub, net.c, net.rxns, net.mets)