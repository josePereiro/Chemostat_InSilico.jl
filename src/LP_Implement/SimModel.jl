mutable struct SimModel
    # net
    net::MetNet              # auxiliar metabolic network
    vatp_idx::Int
    vg_idx::Int
    vl_idx::Int
    obj_idx::Int

    # Sim params
    θvatp::Int               # exactness of vatp discretization dvatp = 10.0^-(θvatp)
    θvg::Int                 # exactness of vg discretization dvatp = 10.0^-(θvg)
    num_min::Float64         # minimal allowed numerical value
    num_max::Float64         # maximal allowed numerical value 
    D::Float64               # Chemostat dilution rate
    
    cg::Float64              # feed G concentration
    Kg::Float64              # g MM constant
    Vg::Float64              # upper g max bound

    cl::Float64              # feed G concentration
    Kl::Float64              # g MM constant
    Vl::Float64 # 0.1        # upper g max bound

    ϵ::Float64               # mutation rate (%)
    niters::Int              # simulation iteration number
    damp::Float64            # numeric damp

    # Chemostat state
    sg::Float64
    sl::Float64
    X::Float64
    Xb::Dict{Float64, Dict{Float64, Float64}}

    function SimModel(;net = ToyModel(), θvatp = 2, θvg = 3, 
            X0 = 0.22, num_min = 1e-30, num_max = 1e30, D = 1e-2, 
            cg = 15.0, sg0 = 15.0, Kg = 0.5, Vg = 0.5, 
            cl = 0.0, sl0 = 0.0, Kl = 0.5, Vl = 0.0, 
            ϵ = 0.0, 
            niters = 200, damp = 0.9,
            sg = Ref{Float64}(sg0), sl = Ref{Float64}(sl0),
            vatp_ider = "vatp",
            vg_ider = "gt",
            vl_ider = "lt",
            obj_ider = "biom"
        )

        vatp_idx = rxnindex(net, vatp_ider)
        vg_idx = rxnindex(net, vg_ider)
        vl_idx = rxnindex(net, vl_ider)
        obj_idx = rxnindex(net, obj_ider)

        Xb = Dict{Float64, Dict{Float64, Float64}}()
        vatp_range, vg_ranges = vatpvg_ranges(net, θvatp, vatp_idx, θvg, vg_idx)
        N = sum(length.(values(vg_ranges)))
        Xi = clamp(X0/N, num_min, num_max)
        for (vatpi, vatp) in enumerate(vatp_range)
            Xb[vatp] = Dict{Float64, Float64}()
            for vg in vg_ranges[vatpi]
                Xb[vatp][vg] = Xi
            end
        end
        X_ts = [X0]
        
        new(
            net, vatp_idx, vg_idx, vl_idx, obj_idx,
            θvatp, θvg, num_min, num_max, D, 
            cg, Kg, Vg, 
            cl, Kl, Vl, 
            ϵ, niters, damp,  
            sg0, sl0, X0, Xb
        )
    end
end

vatpvgN(M::SimModel) = sum(length.(values(M.Xb)))
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