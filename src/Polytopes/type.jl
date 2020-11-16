struct Polytope
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
    # constants
    vgGL::Float64 # vg global lower bound
    vgGU::Float64 # vg global upper bound
    ΔvgG::Float64 # vg global delta
    vatpGL::Float64 # vatp global lower bound
    vatpGU::Float64 # vatp global upper bound
    ΔvatpG::Float64 # vatp global delta


    # Default Constructor
    function Polytope(;Lr = 0.0, Ur = 0.45, Lg = 0.0, Ug = 0.5, Ll = -100.0, Ul = 0.0, 
                Nb = 348.0, Nf = 2, Nr = 18,
                cg = 15, xi = 100.0)
        
        vgGL = Lg
        vgGU = min(Ug, cg / xi)
        vatpGL = Nf * vgGL
        vatpGU = vgGU * Nf * (Nr + 1)

        new(Lr, Ur, Lg, Ug, Ll, Ul, Nb, Nf, Nr, cg, xi, 
            vgGL, vgGU, vgGU - vgGL, vatpGL, vatpGU, vatpGU - vatpGL)
    end
end

is_inpolytope(vatp::Float64, vg::Float64, p::Polytope) =
    vgL(vatp, p) <= vg <= vgU(vatp, p) && vatpL(vg, p) <= vatp <= vatpU(vg, p)
