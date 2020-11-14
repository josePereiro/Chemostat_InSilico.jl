# vg bounds
vg_global_max(p) = min(p.Ug, p.cg / p.xi)
vg_global_max(vatp, p) = vg_global_max(p)
vg_global_min(p) = p.Lg
vg_global_min(vatp, p) = vg_global_min(p)

vg_local_min(vatp, p) = max(vatp/(1 + p.Nr), vatp - p.Ur*p.Nr)/ p.Nf
vg_local_max(vatp, p) = min(vg_global_max(p), vatp / p.Nf)
Δvg(vatp, p) = vg_local_max(vatp, p) - vg_local_min(vatp, p)

##  vatp
vatp_global_min(p) = 0.0
vatp_global_min(vatp, p) = vatp_global_min(p)
vatp_global_max(p) = p.Ur*p.Nr + vg_global_max(p)*p.Nf
vatp_global_max(vatp, p) = vatp_global_max(p)

vatp_local_min(vg, p) = p.Nf * vg
vatp_local_max(vg, p) = min((1 + p.Nr) * p.Nf * vg, p.Ur * p.Nr + p.Nf * vg)
Δvatp(vg, p) = vatp_local_max(vg, p) - vatp_local_min(vg, p)


