rand_vg(p) = (L = vgL(p); U = vgU(p); rand() * (U - L) + L)
rand_vg(vatp, p) = (L = vgL(vatp, p); U = vgU(vatp, p); rand() * (U - L) + L)
rand_vatp(p) = (L = vatpL(p); U = vatpU(p); rand() * (U - L) + L)
rand_vatp(vg, p) = (L = vatpL(vg, p); U = vatpU(vg, p); rand() * (U - L) + L)
