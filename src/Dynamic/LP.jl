const CLP_SOLVER = ClpSolver()
const MAX_SENSE = -1.0
const MIN_SENSE = 1.0

## -------------------------------------------------------------------
function fba(S, b, lb, ub, c; solver = CLP_SOLVER)
    LP = linprog(
        c, # Opt sense vector 
        S, # Stoichiometric matrix
        b, # row lb
        b, # row ub
        lb, # column lb
        ub, # column ub
        solver
    )
    return isempty(LP.sol) ? zeros(length(lb)) : LP.sol
end

function fba(S, b, lb, ub, c, idx, sense; kwargs...)
    bk_c = c[idx]
    c[idx] = sense
    try;
        sol = fba(S, b, lb, ub, c; kwargs...)
        return sol
    finally; c[idx] = bk_c end
end

## -------------------------------------------------------------------
function fva(S, b, lb, ub, c, idx::Integer;  kwargs...)
    bk_c = c[idx]
    c[idx] = MIN_SENSE
    try
        L = fba(S, b, lb, ub, c; kwargs...)[idx]
        c[idx] = MAX_SENSE
        U = fba(S, b, lb, ub, c; kwargs...)[idx]
        return (first(L), first(U))
    finally; c[idx] = bk_c end
end

function fva(S, b, lb, ub, c, idxs; kwargs...)
    L = zeros(length(idxs))
    U = zeros(length(idxs))
    for (i, idx) in idxs |> enumerate
        L[i], U[i] = fva(S, b, lb, ub, c, idx; kwargs...)
    end
    return (L, U)
end