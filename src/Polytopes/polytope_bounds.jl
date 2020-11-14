# vg bounds
vg_global_max(p) = min(p.Ug, p.cg / p.xi)
vg_global_min(p) = p.Lg

vg_local_min(vatp, p) = max(vatp/(1 + p.Nr), vatp - p.Ur*p.Nr)/ p.Nf
vg_local_max(vatp, p) = min(vg_global_max(p), vatp / p.Nf)
Δvg(vatp, p) = vg_local_max(vatp, p) - vg_local_min(vatp, p)

##  vatp
vatp_global_min(p) = p.Nf*vg_global_min(p)
vatp_global_max(p) = vg_global_max(p)*p.Nf*(p.Nr + 1)

vatp_local_min(vg, p) = p.Nf * vg
vatp_local_max(vg, p) = min((1 + p.Nr) * p.Nf * vg, p.Ur * p.Nr + p.Nf * vg)
Δvatp(vg, p) = vatp_local_max(vg, p) - vatp_local_min(vg, p)


