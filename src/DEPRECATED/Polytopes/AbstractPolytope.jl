abstract type AbstractPolytope end

## ------------------------------------------------------------------
is_inpolytope(vatp::Real, vg::Real, p::AbstractPolytope) =
    vgL(vatp, p) <= vg <= vgU(vatp, p) && vatpL(vg, p) <= vatp <= vatpU(vg, p)
