
## ------------------------------------------------------
mutable struct SimD3
    # Space
    V::Space

    # Chemostat
    X::Float64
    sg::Float64
    cg::Float64
    ϵ::Float64
    Δt::Float64
    D::Float64 

    # Simulation
    it::Int
    niters::Int
    dXdt::Float64
    z_av::Float64
    ug_av::Float64
    uo_av::Float64
    cgD_X::Float64
    P::Vector{Float64}

    function SimD3(;
            # Space
            V::Space, 
            # Chemostat
            X::Float64, 
            cg::Float64, 
            sg::Float64, 
            ϵ::Float64, 
            Δt::Float64, 
            D::Float64, 
            # Simulation
            it::Int = 1,
            niters::Int = 5000, 
            dXdt::Float64 = 0.0, 
            z_av::Float64 = 0.0, 
            ug_av::Float64 = 0.0, 
            uo_av::Float64 = 0.0, 
            cgD_X::Float64 = 0.0, 
            P::Vector{Float64} = ones(length(V))
        )

        # Space
        @assert length(P) == length(V)
        normalize!(P)

        new(V, X, sg, cg, ϵ, Δt, D, it, niters, dXdt, z_av, ug_av, uo_av, cgD_X, P)
    end
end

## ------------------------------------------------------
function Base.show(io::IO, S::SimD3)
    println(io, "SimD3")
    println(io, "niters: ", S.niters)
    println(io, S.V)
end
Base.show(S::SimD3) = show(stdout, S)

## ------------------------------------------------------
function run_simD3!(S::SimD3; 
        dobreak::Function,
        tranformation::Function,
        feedback::Function = () -> nothing
    )

    V = S.V 
    P = S.P
    Vz = Vi(V, :z)
    Vug = Vi(V, :ug)
    Vuo = Vi(V, :uo)

    normalize!(P)

    # globals
    cgD_X = 0.0
    dXidt = similar(P)
    loc = similar(P)
    glob = similar(P)
    drain = similar(P)
    vol = length(V)
    
    it0 = (S.it - 1)
    for _ in it0:S.niters
        S.it += 1
        
        # dXidt = loc + glob - drain
        loc .=  (1.0 - S.ϵ) .* Vz .* S.X .* P
        glob .=  S.ϵ * sum(Vz .* S.X .* P) / vol
        drain .=  S.D .* S.X .* P
        dXidt .= (loc .+ glob .- drain) .* S.Δt

        # update P
        P .= (P .+ (dXidt ./ S.X)) 
        normalize!(P)
        S.dXdt = sum(dXidt)
        
        # update X
        S.X += S.dXdt
        
        # P transformation
        S.ug_av = sum(Vug .* P)
        S.cgD_X = S.cg * S.D / S.X
        if (S.ug_av > S.cgD_X)
            tranformation()
            normalize!(P)
        end

        # update exchanges
        S.z_av = sum(Vz .* P)
        S.ug_av = sum(Vug .* P)
        S.uo_av = sum(Vuo .* P)

        # update limiting conc
        S.sg += (- S.ug_av * S.X + (S.cg - S.sg) * S.D) * S.Δt

        # feedback
        feedback()

        # stst
        (dobreak() == true) && break

    end # for it

    return S

end
