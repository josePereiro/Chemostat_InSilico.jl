## ------------------------------------------------------------------
# discretize(val, δ) = (val ÷ δ) * δ
discretize(val, θ::Int) = round(val, digits = θ)
veps(θ::Int) = 10.0^(-θ)

## ------------------------------------------------------------------
struct ΔPolytope <: AbstractPolytope
    # Polytope
    p::Polytope
    # Discretization
    Xs::Dict{Tuple{Float64,Float64}, Float64}
    θvatp::Int
    θvg::Int

    function ΔPolytope(;Lr = 0.0, Ur = 0.45, Lg = 0.0, Ug = 0.5, Ll = -100.0, Ul = 0.0, 
            Nb = 348.0, Nf = 2, Nr = 18,
            D = 0.1, cg = 15, cl = 0.0, X = 10.0,
            θvatp::Int = 2, θvg::Int = 3, 
            X0::Real = 0.0
        ) 
        pol = Polytope(;Lr, Ur, Lg, Ug, Ll, Ul, Nb, Nf, Nr, D, cg, cl, X)
        Xs = initXs(pol, θvatp, θvg, X0)
        new(pol, Xs, θvatp, θvg)
    end
end

## ------------------------------------------------------------------
# ranges

vatp_range(pol::Polytope, θvatp::Int) = vatpL(pol):veps(θvatp):vatpU(pol)
vatp_range(pol::Polytope, θvatp::Int, vg) = vatpL(vg, pol):veps(θvatp):vatpU(vg, pol)
vatp_range(Δp::ΔPolytope) = vatp_range(Δp.p, Δp.θvatp)
vatp_range(Δp::ΔPolytope, vg) = vatp_range(Δp.p, Δp.θvatp, vg)

vatp_nbins(pol::Polytope, θvatp::Int) = length(vatp_range(pol, θvatp))
vatp_nbins(pol::Polytope, θvatp::Int, vg) = length(vatp_range(pol, θvatp, vg))
vatp_nbins(Δp::ΔPolytope) = length(vatp_range(Δp.p, Δp.θvatp))
vatp_nbins(Δp::ΔPolytope, vg) = length(vatp_range(Δp.p, Δp.θvatp, vg))

vg_range(pol::Polytope, θvg::Int) = vgL(pol):veps(θvg):vgU(pol)
vg_range(pol::Polytope, θvg::Int, vatp) = vgL(vatp, pol):veps(θvg):vgU(vatp, pol)
vg_range(Δp::ΔPolytope) = vg_range(Δp.p, Δp.θvg)
vg_range(Δp::ΔPolytope, vatp) = vg_range(Δp.p, Δp.θvg, vatp)

vg_nbins(pol::Polytope, θvg::Int) = length(vg_range(pol, θvg))
vg_nbins(pol::Polytope, θvg::Int, vatp) = length(vg_range(pol, θvg, vatp))
vg_nbins(Δp::ΔPolytope) = length(vg_range(Δp.p, Δp.θvg))
vg_nbins(Δp::ΔPolytope, vatp) = length(vg_range(Δp.p, Δp.θvg, vatp))

nbins(Δp::ΔPolytope) = sum(vg_nbins(Δp, vatp) for vatp in vatp_range(Δp.p, Δp.θvatp))

# ## ------------------------------------------------------------------
# function foreachX(f::Function, Xs, pol, θvatp, θvg)
#     for vatp in vatp_range(pol, θvatp)
#         dvatp = discretize(vatp, θvatp)
#         for vg in vg_range(pol, θvg, vatp)
#             dvg = discretize(vg, θvg) 
#             Xs[(dvatp, dvg)] = f(dvatp, dvg)
#         end
#     end
#     Xs
# end

## ------------------------------------------------------------------
# Setter
setX!(Xs::Dict, vatp::Real, vg::Real, θvatp::Int, θvg::Int, val::Real) = 
    (Xs[(discretize(vatp, θvatp), discretize(vg, θvg))] = val; Xs)
function setX!(f::Function, Xs, pol, θvatp, θvg)
    for vatp in vatp_range(pol, θvatp)
        dvatp = discretize(vatp, θvatp)
        for vg in vg_range(pol, θvg, vatp)
            dvg = discretize(vg, θvg) 
            Xs[(dvatp, dvg)] = f(dvatp, dvg)
        end
    end
    Xs
end
setX!(Δp::ΔPolytope, vatp::Real, vg::Real, val::Real) = setX!(Δp.Xs, vatp, vg, Δp.θvatp, Δp.θvg, val)
setX!(f::Function, Δp::ΔPolytope) = setX!(f, Δp.Xs, Δp.p, Δp.θvatp, Δp.θvg)

## ------------------------------------------------------------------
# Getter
getX(Xs::Dict, vatp::Real, vg::Real, θvatp::Int, θvg::Int, defval::Real) = 
    get(Xs, (discretize(vatp, θvatp), discretize(vg, θvg)), defval)
getX(Xs::Dict, vatp::Real, vg::Real, θvatp::Int, θvg::Int) = 
    Xs[(discretize(vatp, θvatp), discretize(vg, θvg))]
getX(Δp::ΔPolytope, vatp::Real, vg::Real) = getX(Δp.Xs, vatp, vg, Δp.θvatp, Δp.θvg)
getX(Δp::ΔPolytope, vatp::Real, vg::Real, defval::Real) = getX(Δp.Xs, vatp, vg, Δp.θvatp, Δp.θvg, defval)

## ------------------------------------------------------------------
# Init
initXs(pol::Polytope, θvatp::Int, θvg::Int, X0::Real = 0.0) =
    setX!((vatp, vg) -> X0, Dict{Tuple{Float64,Float64}, Float64}(), pol, θvatp, θvg)

## ------------------------------------------------------------------
# merge
function addXs!(Δp1::ΔPolytope, Δp2::ΔPolytope)
    setX!(Δp1) do vatp, vg
        k = (vatp, vg)
        return Δp1.Xs[k] + get(Δp2.Xs, k, 0.0)
    end
    Δp1
end

## ------------------------------------------------------------------
# Independents
# vg
vgL(Δp::ΔPolytope) = vgL(Δp.p)
vgL(vatp, Δp::ΔPolytope) = vgL(vatp, Δp.p)
vgU(Δp::ΔPolytope) = vgU(Δp.p)
vgU(vatp, Δp::ΔPolytope) = vgU(vatp, Δp.p)
Δvg(Δp::ΔPolytope) = Δvg(Δp.p)
Δvg(vatp, Δp::ΔPolytope) = Δvg(vatp, Δp.p)

# Dependents
# vatp
vatpL(Δp::ΔPolytope) = vatpL(Δp.p)
vatpL(vg, Δp::ΔPolytope) = vatpL(vg, Δp.p)
vatpU(vg, Δp::ΔPolytope) = vatpU(vg, Δp.p)
vatpU(Δp::ΔPolytope) = vatpU(Δp.p)
Δvatp(Δp::ΔPolytope) = Δvatp(Δp.p)
Δvatp(vg, Δp::ΔPolytope) = Δvatp(vg, Δp.p)

# vr
vr(vatp, vg, Δp::ΔPolytope) = vr(vatp, vg, Δp.p)

# vl
vl(vatp, vg, Δp::ΔPolytope) = vl(vatp, vg, Δp.p)

# growth
z(vatp, Δp::ΔPolytope) = z(vatp, Δp.p)