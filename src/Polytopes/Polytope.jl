struct Polytope <: AbstractPolytope
    # thermodinamic bounds
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
    D::Float64
    cg::Float64
    cl::Float64
    X::Float64
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
                D = 0.1, cg = 15, cl = 0.0, X = 10.0
            )
        
        vgGL = Lg
        vgGU = vgU(Ug, cg, X, D)
        vatpGL = Nf * vgGL
        vatpGU = vgGU * Nf * (Nr + 1)

        new( 
            Lr, Ur, Lg, Ug, Ll, Ul, 
            Nb, Nf, Nr, 
            D, cg, cl, X, 
            vgGL, vgGU, vgGU - vgGL, vatpGL, vatpGU, vatpGU - vatpGL
        )
    end
end


# TODO: Complete all fluxes bounds
# TODO: Make possible to eat lactate

# Independents
# vg
vgL(p::Polytope) = p.vgGL
vgL(vatp, p::Polytope) = max(vatp/(1 + p.Nr), vatp - p.Ur * p.Nr)/ p.Nf
vgU(Ug, cg, X, D) = min(Ug, (cg * D) / X)
vgU(p::Polytope) = p.vgGU
vgU(vatp, p::Polytope) = min(vgU(p), vatp / p.Nf)
Δvg(p::Polytope) = p.ΔvgG
Δvg(vatp, p::Polytope) = vgU(vatp, p) - vgL(vatp, p)

# Dependents
# vatp
vatpL(p::Polytope) = p.vatpGL
vatpL(vg, p::Polytope) = p.Nf * vg
vatpU(vg, p::Polytope) = min((1 + p.Nr) * p.Nf * vg, p.Ur * p.Nr + p.Nf * vg)
vatpU(p::Polytope) = p.vatpGU
Δvatp(p::Polytope) = p.ΔvatpG
Δvatp(vg, p::Polytope) = vatpU(vg, p) - vatpL(vg, p)

# vr
vr(vatp, vg, p::Polytope) = (vatp - p.Nf * vg)/ p.Nr

# vl
vlU(Ul, cl, X, D) = min(Ul, (cl * D) / X)
vl(vatp, vg, p::Polytope) = vr(vatp, vg, p) - p.Nf * vg

# growth
z(vatp, p::Polytope) = vatp / p.Nb

