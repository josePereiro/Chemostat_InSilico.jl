## ----------------------------------------------------------------------------
function _fba_max_vatp_min_G(M; δ, LP_cache, verbose)

    vatp_range, vg_ranges = Dyn.vatpvg_ranges(M)
    max_vatp = -Inf
    min_vg = Inf
    # (max_vatp, min_vg) will maximize the yield
    for (vatpi, vatp) in vatp_range |> enumerate
        vg_range = vg_ranges[vatpi]
        isempty(vg_range) && continue
        if vatp > max_vatp
            max_vatp = vatp
            min_vg = minimum(vg_range)
        end
    end
    @assert !isinf(max_vatp)

    fbaf(vatp, vg) = (vatp == max_vatp && vg == min_vg) ? 1.0 : 0.0
    FBAMs = Dyn.get_marginals(fbaf, M; δ, LP_cache, verbose)

    return FBAMs
end

function _fba_max_vatp_max_G(M; δ, LP_cache, verbose)

    vatp_range, vg_ranges = Dyn.vatpvg_ranges(M)
    max_vatp = -Inf
    max_vg = -Inf
    # (max_vatp, min_vg) will maximize the yield
    for (vatpi, vatp) in vatp_range |> enumerate
        vg_range = vg_ranges[vatpi]
        isempty(vg_range) && continue
        if vatp > max_vatp
            max_vatp = vatp
            max_vg = maximum(vg_range)
        end
    end
    @assert !isinf(max_vatp)

    fbaf(vatp, vg) = (vatp == max_vatp && vg == max_vg) ? 1.0 : 0.0
    FBAMs = Dyn.get_marginals(fbaf, M; δ, LP_cache, verbose)

    return FBAMs
end

## ----------------------------------------------------------------------------
function run_FBA!(M, FBAmode; LP_cache, δ, δμ, biom_avPX, verbose = true)

    # Z FIXXED
    ismode = FBAmode in [FBA_Z_FIXXED_G_OPEN, FBA_Z_FIXXED_G_BOUNDED, FBA_Z_FIXXED_MAX_G, FBA_Z_FIXXED_MIN_G]
    ismode && let
        # Fix biomass to observable
        net = M.net
        net.ub[M.obj_idx] = biom_avPX * (1.0 + δμ)
        net.lb[M.obj_idx] = biom_avPX * (1.0 - δμ)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    # G BOUNDING
    ismode = FBAmode in [FBA_Z_OPEN_G_BOUNDED, FBA_Z_FIXXED_G_BOUNDED, FBA_Z_FIXXED_MAX_G, FBA_Z_FIXXED_MIN_G]
    ismode && let
        # Fix av_ug
        net = M.net
        net.ub[M.vg_idx] = min(M.Vg, M.cg * M.D/ M.X)
        net.ub[M.vl_idx] = min(M.Vl, M.cl * M.D/ M.X)
        L, U = Dyn.fva(net)
        net.lb .= L; net.ub .= U
    end

    # FBA
    ismode = FBAmode in [FBA_Z_FIXXED_MAX_G]
    if ismode
        FBAMs = _fba_max_vatp_max_G(M; δ, LP_cache, verbose)
    else
        FBAMs = _fba_max_vatp_min_G(M; δ, LP_cache, verbose)
    end
    return FBAMs
end
