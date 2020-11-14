struct Cell
    p::SimParams
    vatp::Float64
    vg::Float64
end

vatp(c::Cell) = c.vatp
vg(c::Cell) = c.vg
# TODO: finish other dependen fluxes

is_inpolytope(c::Cell) = is_inpolytope(c.vatp, c.vg, c.p)
