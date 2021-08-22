# Store network information
# TODO: echenolize model

struct MetNet

    ## LP (original) 
    S::Matrix{Float64}
    b::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    c::Vector{Float64}
    rxns::Vector{String}
    mets::Vector{String}

end

MetNet(;S, b, lb, ub, c, rxns, mets) = MetNet(S, b, lb, ub, c, rxns, mets)

## -------------------------------------------------------------------
# Utils
_IDER_TYPE = Union{String, Symbol, Int}
_IDER_TYPE_STR = Union{String, Symbol}
_IDER_TYPE_INT = Int

rxnindex(net::MetNet, rxn::_IDER_TYPE_STR) = findfirst(isequal(string(rxn)), net.rxns)
rxnindex(net::MetNet, rxn::_IDER_TYPE_INT) = findfirst(isequal(rxn), eachindex(net.rxns))
rxnindex(net::MetNet, rxns::Vector) = rxnindex.([net], rxns)
rxnindex(net::MetNet, rxns) = rxnindex(net, collect(rxns))

metindex(net::MetNet, met::_IDER_TYPE_STR) = findfirst(isequal(string(met)), net.mets)
metindex(net::MetNet, met::_IDER_TYPE_INT) = findfirst(isequal(met), eachindex(net.mets))
metindex(net::MetNet, mets::Vector) = metindex.([net], mets)
metindex(net::MetNet, mets) = metindex(net, collect(mets))

rxns(net::MetNet, rxns) = net.rxns[rxnindex(net, rxns)]

## -------------------------------------------------------------------
function reducebox!(net::MetNet)
    lb, ub = fva(net)
    net.lb .= lb
    net.ub .= ub
    return net
end
_reducebox!(net, reduce) = reduce ? reducebox!(net) : net

function fix!(net::MetNet, rxns, val; reducebox = false)
    reducebox && reducebox!(net)
    idxs = rxnindex(net, rxns)
    idxs = (idxs isa Vector) ? idxs : [idxs]
    net.lb[idxs] .= val
    net.ub[idxs] .= val
    net
end

function fixxing(f::Function, net::MetNet, rxns, val; tol = 0.0)
    
    idxs = rxnindex(net, rxns)
    bk_lb = net.lb[idxs]
    bk_ub = net.ub[idxs]
    net.lb[idxs] = (val .- tol)
    net.ub[idxs] = (val .+ tol)
    try; return f()
    finally; 
        net.lb[idxs] = bk_lb
        net.ub[idxs] = bk_ub
    end
end

rxnid(net::MetNet, rxn) = net.rxn[rxnindex(net, rxn)]
lb(net::MetNet, rxn) = net.lb[rxnindex(net, rxn)]
lb!(net::MetNet, rxn, val; reducebox = false) = 
    (net.lb[rxnindex(net, rxn)] = val; _reducebox!(net, reducebox))
ub(net::MetNet, rxn) = net.ub[rxnindex(net, rxn)]
ub!(net::MetNet, rxn, val; reducebox = false) = 
    (net.ub[rxnindex(net, rxn)] = val; _reducebox!(net, reducebox))
bounds!(net::MetNet, rxn, lb, ub; reducebox = false) = 
    (idx = rxnindex(net, rxn); net.lb[idx] = lb; net.ub[idx] = ub; _reducebox!(net, reducebox))
bounds(net::MetNet, rxn) = 
    (idx = rxnindex(net, rxn); (net.lb[idx],  net.ub[idx]))

## -------------------------------------------------------------------
_get_elstr(v, i) = i > lastindex(v) ? "" : string(v[i])
function pretty(;kwargs...)
    kwargs = collect(kwargs)
    ks = first.(kwargs)
    vs = last.(kwargs)
    # maximum length
    padl = map(zip(ks, vs)) do (k, v)
        maxvl = maximum(length.(string.(v)))
        max(maxvl, length(string(k)))
    end
    
    tab = "   "
    println(join(rpad.(first.(kwargs), padl), tab))
    println(join( ["-"].^padl , tab))
    maxl = maximum(length.(last.(kwargs)))
    for i in 1:maxl
        println(join(rpad.(_get_elstr.(vs, i), padl), tab))
    end
    
end

## -------------------------------------------------------------------
_round(n) = round(n; sigdigits = 3)
function prettyrxn(net::MetNet, v::Vector, rxns = eachindex(net.rxns))
    rxns = rxnindex(net, rxns)
    pretty(;
        n = collect(eachindex(net.rxns))[rxns], 
        rxn = net.rxns[rxns], 
        v = _round.(v)[rxns], 
        lb = _round.(net.lb)[rxns], 
        ub = _round.(net.ub)[rxns], 
        eq = rxn_str.([net], net.rxns)[rxns]
    )
end

## -------------------------------------------------------------------
# LP
fba(net::MetNet; kwargs...) = fba(net.S, net.b, net.lb, net.ub, net.c; kwargs...)
fba(net::MetNet, obj, sense = MAX_SENSE; kwargs...) = 
    fba(net.S, net.b, net.lb, net.ub, net.c, rxnindex(net, obj), sense; kwargs...)
fva(net::MetNet, rxns = eachindex(net.rxns); kwargs...) = 
    fva(net.S, net.b, net.lb, net.ub, net.c, rxnindex(net, rxns); kwargs...)
Δv(net, rxns = eachindex(net.rxns)) = 
    ((lb, ub) = fva(net, rxnindex(net, rxns)); ub - lb)
U(net, rxns = eachindex(net.rxns)) = 
    ((lb, ub) = fva(net, rxnindex(net, rxns)); ub)
L(net, rxns = eachindex(net.rxns)) = 
    ((lb, ub) = fva(net, rxnindex(net, rxns)); lb)

## -------------------------------------------------------------------
Base.hash(net::MetNet) = hash((:MetNet, net.S, net.b, net.lb, net.ub))
Base.size(net::MetNet) = size(net.S)

## -------------------------------------------------------------------
function rxn_str(net::MetNet, rxn)
    ri = rxnindex(net, rxn)
    rS = net.S[:, ri]
    
    reacts = String[]
    prods = String[]
    for (mi, ms) in enumerate(rS)
        mid = net.mets[mi]
        ms = round(ms; sigdigits = 3)
        (ms < 0.0) && push!(reacts, string("(", ms, ")", mid))
        (ms > 0.0) && push!(prods, string("(", ms, ")", mid))
    end

    reacts_str = join(reacts, " + ")
    prods_str = join(prods, " + ")

    lb, ub = net.lb[ri], net.ub[ri]
    arrow = (lb == ub == 0) ? 
        ">-<" : (lb < 0.0 && ub > 0.0) ? 
            "<-->" : (lb >= 0.0 && ub > 0.0) ? 
                "-->" : (lb < 0.0 && ub <= 0.0) ? 
                    "<--" : "error"
    
    
    string(reacts_str, " ", arrow, " ", prods_str)
    
end

## -------------------------------------------------------------------
function Base.show(io::IO, net::MetNet)
    
    println(io, "MetNet: ", size(net))

    M, N = size(net)

    ri_strs = String["ri"; string.(1:N)]
    rxnid_strs = String["rxnid"]
    bounds_strs = String["bounds"]
    rxni_strs = String["rxn str"]
    for ri in 1:N
        lb, ub = net.lb[ri], net.ub[ri]
        lb = round(lb; sigdigits = 3)
        ub = round(ub; sigdigits = 3)
        bounds_str = string("[", lb, ", " ,ub, "]")
        push!(bounds_strs, bounds_str)
        
        rid = net.rxns[ri]
        rxnid_str = string(rid)
        push!(rxnid_strs, rxnid_str)
        
        rxni_str = rxn_str(net, ri)
        push!(rxni_strs, rxni_str)

    end

    ri_pad = maximum(length.(ri_strs))
    rxnid_pad = maximum(length.(rxnid_strs))
    rxni_pad = maximum(length.(rxni_strs))
    bounds_pad = maximum(length.(bounds_strs))
    tab = "   "
    for ri in eachindex(rxnid_strs)
        println(io, 
            join([
                rpad(ri_strs[ri], ri_pad),  
                rpad(rxnid_strs[ri], rxnid_pad),  
                rpad(bounds_strs[ri], bounds_pad), 
                rpad(rxni_strs[ri], rxni_pad)
            ], tab)
        )
    end
    nothing
end

Base.show(net::MetNet) = show(stdout, net)

## -------------------------------------------------------------------
function is_feasible(net::MetNet, freeis, v; tol::Vector = Δv(net, freeis) .* 0.01)
    fullv = fixxing(net, freeis, v; tol) do
        fba(net)
    end
    any(fullv .!= 0.0)
end

## -------------------------------------------------------------------
function BoxGrid(net::MetNet, frees::Vector, δs::Vector)
    idxs = rxnindex.([net], frees)
    dimdat = []
    for (idx, free, δ) in zip(idxs, frees, δs)
        isnothing(idx) && error(free, " not found")
        id = net.rxns[idx]
        lb, ub = net.lb[idx], net.ub[idx]
        push!(dimdat, (Symbol(id), lb, ub, δ))
    end
    BoxGrid(dimdat)
end
BoxGrid(net::MetNet, frees::Vector, δ::Int) = BoxGrid(net, frees, fill(δ, length(frees)))
BoxGrid(net::MetNet, δ::Int) = BoxGrid(net, net.rxns, fill(δ, length(net.rxns)))
