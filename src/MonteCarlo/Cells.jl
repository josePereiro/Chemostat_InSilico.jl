struct Cell
    p::Polytope
    vatp::Float64
    vg::Float64
end

vatp(c::Cell) = c.vatp
vg(c::Cell) = c.vg
# TODO: finish other dependen fluxes

Polytopes.is_inpolytope(c::Cell) = Polytopes.is_inpolytope(c.vatp, c.vg, c.p)