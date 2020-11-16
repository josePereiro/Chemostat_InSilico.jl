# vg bounds
vgL(p::Polytope) = p.vgGL
vgL(vatp, p::Polytope) = max(vatp/(1 + p.Nr), vatp - p.Ur * p.Nr)/ p.Nf
vgU(p::Polytope) = p.vgGU
vgU(vatp, p::Polytope) = min(vgU(p), vatp / p.Nf)
Δvg(p::Polytope) = p.ΔvgG
Δvg(vatp, p::Polytope) = vgU(vatp, p) - vgL(vatp, p)

##  vatp
vatpL(p::Polytope) = p.vatpGL
vatpL(vg, p::Polytope) = p.Nf * vg
vatpU(vg, p::Polytope) = min((1 + p.Nr) * p.Nf * vg, p.Ur * p.Nr + p.Nf * vg)
vatpU(p::Polytope) = p.vatpGU
Δvatp(p::Polytope) = p.ΔvatpG
Δvatp(vg, p::Polytope) = vatpU(vg, p) - vatpL(vg, p)

