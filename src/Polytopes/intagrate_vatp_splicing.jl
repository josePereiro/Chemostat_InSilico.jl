function integrate_vatp_slicing(f::Function, p, 
    vatp0 = vatp_global_min(p), vatp1 = vatp_global_max(p); n = Int(1e3))
    
    @assert vatp0 < vatp1 string("vatp0(", vatp0, ") < vatp1(", vatp1, ")")
    rvatp = range(vatp0, vatp1; length = n)
    Δvatp = step(rvatp)
    res = 0.0
    for vatp in rvatp[1:end-1]
        res += f(vatp) * Δvg(vatp, p) * Δvatp
    end
    res
end
