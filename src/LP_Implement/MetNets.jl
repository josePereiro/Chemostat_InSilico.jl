# Store network information
# TODO: echenolize model

const ABS_MAX_BOUND = 1000.0 # abs bound
const BIOMASS_IDER = "biom"
const RXNS = [ "gt"  , "ferm" , "resp" , "ldh" ,  "lt" , BIOMASS_IDER , "atpm" , "vatp" ]
const METS = ["G", "E", "P", "L", "AUX1"]
const Y = 348 

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
function ToyModel()
    net = Dict()

    # Network parameters from 
    # Fernandez-de-Cossio-Diaz, Jorge, Roberto Mulet, and Alexei Vazquez. (2019) https://doi.org/10.1038/s41598-019-45882-w.
    Nf = 2.0
    Nr = 38.0
    y = -Y        # mmol/gDW
    Vr = 0.45     # mmol/gDW h
    ATPM = 0.0 # 1.0625 # mmol/gDW h
    Vg = 0.5      # mmol/gDW h

    net[:S] = 
    # rxns: gt    ferm  resp  ldh   lt   biom    atpm  vatp    # mets
    [       1.0  -1.0   0.0   0.0   0.0   0.0    0.0   0.0  ;  #  G
            0.0    Nf    Nr   0.0   0.0    y    -1.0   0.0  ;  #  E
            0.0    Nf  -1.0  -1.0   0.0   0.0    0.0   0.0  ;  #  P
            0.0   0.0   0.0   1.0   1.0   0.0    0.0   0.0  ;  #  L
            0.0    Nf    Nr   0.0   0.0   0.0    0.0  -1.0  ;  #  AUX1
    ]
    
    net[:mets] = deepcopy(METS)
    net[:b] =    [0.0, 0.0, 0.0, 0.0, 0.0] # const exchanges
    
    AB = ABS_MAX_BOUND
    net[:rxns] = deepcopy(RXNS)
    #            [ "gt"  , "ferm" , "resp" , "ldh" ,  "lt" , BIOMASS_IDER , "atpm" , "vatp" ]
    net[:lb]   = [ 0.0   ,  0.0   ,  0.0   ,  0.0  ,  -AB  ,     0.0      ,  ATPM  ,  0.0   ];
    net[:ub]   = [  Vg   ,  AB    ,   Vr   ,  AB   ,   0.0 ,      AB      ,   AB   ,  AB    ];
    net[:c]    = [ 0.0   ,  0.0   ,  0.0   ,  0.0  ,   0.0 ,  MAX_SENSE   ,   0.0  ,  0.0   ];
    return MetNet(;net...)
end

## -------------------------------------------------------------------
# LP
fba(net::MetNet; solver = CLP_SOLVER) = fba(net.S, net.b, net.lb, net.ub, net.c; solver)
fba(net::MetNet, objidx, sense = MAX_SENSE; solver = CLP_SOLVER) = 
    fba(net.S, net.b, net.lb, net.ub, net.c, objidx, sense; solver)
fva(net::MetNet, idx::Integer; solver = CLP_SOLVER) = 
    fva(net.S, net.b, net.lb, net.ub, net.c, idx; solver)
fva(net::MetNet, idxs = eachindex(net.rxns); solver = CLP_SOLVER) = 
    fva(net.S, net.b, net.lb, net.ub, net.c, idxs; solver)
Δv(net, idx) = begin lb, ub = fva(net, idx); ub - lb end
U(net, idx) = begin lb, ub = fva(net, idx); ub end
L(net, idx) = begin lb, ub = fva(net, idx); lb end

## -------------------------------------------------------------------
# Utils
rxnindex(net::MetNet, id) = findfirst(isequal(id), net.rxns)
metindex(net::MetNet, id) = findfirst(isequal(id), net.rxns)
fix!(net, idx, val) = net.lb[idx] = net.ub[idx] = val
function fixxing(f, net, idx, val)
    bk_lb, bk_ub = net.lb[idx], net.ub[idx]
    net.lb[idx] = net.ub[idx] = val
    val = f()
    net.lb[idx], net.ub[idx] = bk_lb, bk_ub
    return val
end

## -------------------------------------------------------------------
Base.hash(net::MetNet) = hash((net.S, net.b, net.lb, net.ub))