mutable struct SimModel
    # net
    net::MetNet              # metabolic network
    vatp_idx::Int
    vg_idx::Int
    vl_idx::Int
    obj_idx::Int

    # Sim params
    θvatp::Int               # exactness of vatp discretization dvatp = 10.0^-(θvatp)
    θvg::Int                 # exactness of vg discretization dvatp = 10.0^-(θvg)
    
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
    function SimModel(;net = ToyModel(), θvatp = 2, θvg = 3, Δt = 0.1,
            X0 = 0.22, D0 = 1e-2,                                          
            cg = 15.0, sg0 = 15.0, Kg = 0.5, Vg = 0.5,
            cl = 0.0, sl0 = 0.0, Kl = 0.5, Vl = 0.1,
            σ = 0.01, τ = 0.0022, ϵ = 0.0,                                            
            num_min = -1e30, num_max = 1e30, 
            niters = 20000, 
            vatp_ider = "vatp", vg_ider = "gt", vl_ider = "lt", obj_ider = "biom"
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
        vatp_range, vg_ranges = vatpvg_ranges(net, θvatp, vatp_idx, θvg, vg_idx)
        i_vatp_range = enumerate(vatp_range)
        N = sum(length.(values(vg_ranges)))
        Xi = clamp(X0/N, 0.0, num_max)
        complete_board!(Xb, i_vatp_range, vg_ranges, Xi)
        X_ts = [X0]
        
        new(
            net, vatp_idx, vg_idx, vl_idx, obj_idx, 
            θvatp, θvg, 
            cg, Kg, Vg, 
            cl, Kl, Vl, 
            ϵ, σ, τ,
            Δt, niters,
            D0, sg0, sl0, X0, Xb
        )
    end
end

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

function ploy_box(M::SimModel)
    vatp_range, vg_ranges = vatpvg_ranges(M)
    vatpL, vatpU = first(vatp_range), last(vatp_range)
    vgL, vgU = minimum(first.(vg_ranges)), maximum(last.(vg_ranges))
    vatpL, vatpU, vgL, vgU
end