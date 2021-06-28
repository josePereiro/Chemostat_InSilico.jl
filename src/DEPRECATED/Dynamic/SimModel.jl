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
    Vg::Float64              # (mmol/ gDW h) upper glc max bound

    cl::Float64              # (mM) feed lac concentration
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
            cg = 15.0, sg0 = 15.0, Vg = 0.5,
            cl = 0.0, sl0 = 0.0, Vl = 0.0, # 0.1
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
        net.ub[vg_idx] = max(net.lb[vg_idx], Vg)
        net.ub[vl_idx] = max(net.lb[vl_idx], Vl)

        # board
        Xb = Dict{Float64, Dict{Float64, Float64}}()
        if fill_Xb
            vatp_range, vg_ranges = vatpvg_ranges(net, δvatp, vatp_idx, δvg, vg_idx)
            i_vatp_range = enumerate(vatp_range)
            N = sum(length.(values(vg_ranges)))
            Xi = clamp(X0/N, 0.0, num_max)
            complete_board!(Xb, i_vatp_range, vg_ranges, Xi)
        end    
        
        
        M = new(
            net, vatp_idx, vg_idx, vl_idx, obj_idx, 
            δvatp, δvg, 
            cg, Vg, 
            cl, Vl, 
            ϵ, σ, τ,
            Δt, niters,
            D0, sg0, sl0, X0, Xb
        )
        check_params(M)
        return M
    end
end

Base.merge!(M1, M2) = foreach((f) -> setfield!(M1, f, getfield(M2, f)), fieldnames(SimModel))
function Base.hash(M::SimModel)
    fs = filter(!isequal(:Xb), fieldnames(typeof(M)))
    hash(string(getproperty.([M], fs)))
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


const CG_MAX_BOUNDS = (0.0, 15.0)
const VG_MAX_BOUNDS = (0.5, 0.5)
const CL_MAX_BOUNDS = (0.0, 0.0)
const VL_MAX_BOUNDS = (0.0, 0.0)
const D_MAX_BOUNDS = (0.0, 0.05)
const X_MAX_BOUNDS = (0.01, 5.0)
const REF_METNET = ToyModel()

function check_params(M::SimModel)

    # field
    for (field, (lb, ub)) in [
            (:cg, CG_MAX_BOUNDS), (:Vg, VG_MAX_BOUNDS),
            (:cl, CL_MAX_BOUNDS), (:Vl, VL_MAX_BOUNDS),
            (:D, D_MAX_BOUNDS), (:X, X_MAX_BOUNDS),
        ]
        val = getproperty(M, field)
        !(lb <= val <= ub) && error(
            string(field, "=", val, ", It must be between ", (lb, ub))
        )
    end

    # net
    for rxni in eachindex(REF_METNET.rxns)
        (REF_METNET.lb[rxni] > M.net.lb[rxni] || REF_METNET.ub[rxni] < M.net.ub[rxni]) && error(
            string(M.net.rxns[rxni], " bounds=", (M.net.lb[rxni], M.net.ub[rxni]), 
            ", It must be between ", (REF_METNET.lb[rxni] , REF_METNET.ub[rxni]))
        )
    end
end

drop_Xb!(M::SimModel) = (empty!(M.Xb); M)