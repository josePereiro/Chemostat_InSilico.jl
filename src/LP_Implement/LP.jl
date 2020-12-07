const CLP_SOLVER = ClpSolver()
const MAX_SENSE = -1.0
const MIN_SENSE = 1.0

## -------------------------------------------------------------------
function fba(S, b, lb, ub, c; solver = CLP_SOLVER)

    sol = linprog(
        c, # Opt sense vector 
        S, # Stoichiometric matrix
        b, # row lb
        b, # row ub
        lb, # column lb
        ub, # column ub
        solver
    )
    return sol.sol
end

function fba(S, b, lb, ub, c, idx, sense; solver = CLP_SOLVER)
    bk_c = c[idx]
    c[idx] = sense
    sol = fba(S, b, lb, ub, c; solver)
    c[idx] = bk_c
    return sol
end

## -------------------------------------------------------------------
function fva(S, b, lb, ub, c, idx::Integer; solver = CLP_SOLVER)
    bk_c = c[idx]
    c[idx] = MIN_SENSE
    L = fba(S, b, lb, ub, c; solver)[idx]
    c[idx] = MAX_SENSE
    U = fba(S, b, lb, ub, c; solver)[idx]
    c[idx] = bk_c
    return (L, U)
end

function fva(S, b, lb, ub, c, idxs; solver = CLP_SOLVER)
    L = zeros(length(idxs))
    U = zeros(length(idxs))
    for (i, idx) in idxs |> enumerate
        L[i], U[i] = fva(S, b, lb, ub, c, idx; solver)
    end
    return (L, U)
end