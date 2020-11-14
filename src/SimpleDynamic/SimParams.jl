# Constants
struct SimParams
    # bounds
    Lr::Float64
    Ur::Float64
    Lg::Float64
    Ug::Float64
    Ll::Float64
    Ul::Float64
    # stoi coe
    Nb::Float64
    Nf::Float64
    Nr::Float64
    # chemostat
    cg::Float64
    xi::Float64
    SimParams(;Lr = 0.0, Ur = 0.45, Lg = 0.0, Ug = 0.5, Ll = -100.0, Ul = 0.0, 
                Nb = 348, Nf = 2, Nr = 18,
                cg = 15, xi = 1) = 
            new(Lr, Ur, Lg, Ug, Ll, Ul, Nb, Nf, Nr, cg, xi)
end

