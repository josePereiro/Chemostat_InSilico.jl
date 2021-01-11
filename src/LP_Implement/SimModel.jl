mutable struct SimModel
    # net
    net::MetNet              # metabolic network

    # important indexes
    vatp_idx::Int
    vg_idx::Int
    vl_idx::Int
    obj_idx::Int

    # Sim params
    δvatp::Int               # exactness of vatp discretization dvatp = 10.0^-(δvatp)
    δvg::Int                 # exactness of vg discretization dvatp = 10.0^-(δvg)
    
    cg::Float64              # (mM) feed glc concentration
    Kg::Float64              # (mM) glc MM constant
    Vg::Float64              # (mmol/ gDW h) upper glc max bound

    cl::Float64              # (mM) feed lac concentration
    Kl::Float64              # (mM) lac MM constant
    Vl::Float64              # (mmol/ gDW h) upper lac max bound

    ϵ::Float64               # mutation rate (%)
    σ::Float64               # (1/ h) dead fraction
    τ::Float64               # (1/ h mM) toxicity coefficient 
    
    Δt::Float64              # (h) simulation time step
    niters::Int              # simulation iteration number

    # Chemostat state
    D::Float64               # (1/ h) Chemostat dilution rate
    sg::Float64              # (mM) medium glc concentration
    sl::Float64              # (mM) medium lac concentration
    X::Float64               # (gDW/ L) cell concentration
    Xb::Dict{Float64, Dict{Float64, Float64}}

    # Default values from
    # Fernandez-de-Cossio-Diaz, Jorge, Roberto Mulet, and Alexei Vazquez. (2019) https://doi.org/10.1038/s41598-019-45882-w.
    function SimModel(;net = ToyModel(), δvatp = 2, δvg = 3, Δt = 0.1,
            X0 = 1.5, D0 = 1e-2,                                          
            cg = 15.0, sg0 = 15.0, Kg = 0.5, Vg = 0.5,
            cl = 0.0, sl0 = 0.0, Kl = 0.5, Vl = 0.1,
            σ = 0.01, τ = 0.0022, ϵ = 0.0,                                            
            num_min = -1e30, num_max = 1e30, 
            niters = 20000, 
            vatp_ider = "vatp", vg_ider = "gt", vl_ider = "lt", obj_ider = "biom",
            fill_Xb = true
        )

        # indexes
        vatp_idx = rxnindex(net, vatp_ider)
        vg_idx = rxnindex(net, vg_ider)
        vl_idx = rxnindex(net, vl_ider)
        obj_idx = rxnindex(net, obj_ider)

        # update net
        net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * sg0) / (Kg + sg0))
        net.ub[vl_idx] = max(net.lb[vl_idx], (Vl * sl0) / (Kl + sl0))

        # board
        Xb = Dict{Float64, Dict{Float64, Float64}}()
        if fill_Xb
            vatp_range, vg_ranges = vatpvg_ranges(net, δvatp, vatp_idx, δvg, vg_idx)
            i_vatp_range = enumerate(vatp_range)
            N = sum(length.(values(vg_ranges)))
            Xi = clamp(X0/N, 0.0, num_max)
            complete_board!(Xb, i_vatp_range, vg_ranges, Xi)
        end    
        X_ts = [X0]
        
        new(
            net, vatp_idx, vg_idx, vl_idx, obj_idx, 
            δvatp, δvg, 
            cg, Kg, Vg, 
            cl, Kl, Vl, 
            ϵ, σ, τ,
            Δt, niters,
            D0, sg0, sl0, X0, Xb
        )
    end
end

Base.merge!(M1, M2) = foreach((f) -> setfield!(M1, f, getfield(M2, f)), fieldnames(SimModel))

vatpvgN(Xb) = sum(length.(values(Xb)))
vatpvgN(M::SimModel) = vatpvgN(M.Xb)

function lXgamma(M::SimModel)
    mX, MX = Inf, -Inf
    for (vatp, Xvatp) in M.Xb
        for (vg, X) in Xvatp
            (X > MX) && (MX = X)
            (X < mX) && (mX = X)
        end
    end
    return mX, MX
end

function pol_box(M::SimModel)
    vatp_range, vg_ranges = vatpvg_ranges(M)
    vatpL, vatpU = first(vatp_range), last(vatp_range)
    vgL, vgU = minimum(first.(vg_ranges)), maximum(last.(vg_ranges))
    (;vatpL, vatpU, vgL, vgU)
end

function plot_box_grid(M::SimModel)
    vatpL, vatpU, vgL, vgU = pol_box(M)
    vatp_range = vrange(vatpL, vatpU, M.δvatp)
    vg_range = vrange(vgL, vgU, M.δvg)
    vatp_range, vg_range
end
