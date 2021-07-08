
## ------------------------------------------------------
mutable struct Sim2D
    # Space
    V::NamedTuple

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
    ug_av::Float64
    z_av::Float64
    P::Vector{Float64}

    function Sim2D(;
            # Space
            Vcell::Space, 
            # Chemostat
            X::Float64, 
            sg::Float64, 
            cg::Float64, 
            ϵ::Float64, 
            Δt::Float64, 
            D::Float64, 
            # Simulation
            it::Int = 1, 
            niters::Int = 5000, 
            dXdt::Float64 = 0.0, 
            ug_av::Float64 = 0.0, 
            z_av::Float64 = 0.0, 
            P::Vector{Float64} = ones(length(Vcell))
        )

        # Space
        vol = length(Vcell)

        z = Vi(Vcell, :z)
        Uz = maximum(z)
        Lz = minimum(z)

        ug = Vi(Vcell, :ug)
        Uug = maximum(ug)
        Lug = minimum(ug)
        
        V = (;vol, z, Uz, Lz, ug, Uug, Lug)
        
        @assert length(P) == vol
        normalize!(P)

        new(V, X, sg, cg, ϵ, Δt, D, it, niters, dXdt, ug_av, z_av, P)
    end

end

## ------------------------------------------------------
function run_simD2!(S; 
        dobreak::Function,
        tranformation::Function,
        feedback::Function = () -> nothing
    )

    @extract S: V P

    # globals
    cgD_X = 0.0
    dXidt = similar(P)
    
    it0 = S.it - 1
    for _ in it0:S.niters
        S.it += 1
            
        # dXidt
        loc =  (1.0 - S.ϵ) .* V.z .* S.X .* P
        glo =  S.ϵ * sum(V.z .* S.X .* P) / V.vol
        drain =  S.D .* S.X .* P
        dXidt .= (loc .+ glo .- drain) .* S.Δt

        # update
        P .= (P .+ (dXidt ./ S.X)) 
        normalize!(P)
        S.dXdt = sum(dXidt)
        S.X += S.dXdt
        
        S.ug_av = sum(V.ug .* P)
        S.z_av = sum(V.z .* P)

        S.sg += (- S.ug_av * S.X + (S.cg - S.sg) * S.D) * S.Δt
        
        # P transformation
        cgD_X = S.cg * S.D / S.X
        if (S.ug_av > cgD_X)
            tranformation()
            normalize!(P)
            S.ug_av = sum(V.ug .* P)
        end

        # feedback
        feedback()

        # stst
        (dobreak() == true) && break

    end # for it

    return S

end
