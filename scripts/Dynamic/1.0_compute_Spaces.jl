using ProjAssistant
@quickactivate

@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic

    using Base.Threads
    using Plots
    import GR
    !isinteractive() && GR.inline("png")
end


## ------------------------------------------------------
# common globals
let
    sglob(Dyn;
        δ = 400
    )
end

## ------------------------------------------------------
@time let
    # input
    net2D = Dyn.ToyModel2D()
    
    freeids = [:z, :ug]
    @lglob Dyn δ
    
    @info("Building Vcell2D", nthreads())
    @show freeids δ
    filter(net, v) = all(v .!= 0.0) # avoid artifacts at 0.0
    Vcell2D = Dyn.Space(net2D, freeids, δ; filter)
    @sglob Dyn Vcell2D net2D

    ps = Dyn.pol_plots(Vcell2D; bins = 10000)
    sfig(Dyn, ps, "Vcell2D.png")
end
 
## ------------------------------------------------------
@time let
    # input
    net3D = Dyn.ToyModel3D()
    
    freeids = [:z, :ug, :uo]
    @lglob Dyn δ
    
    @info("Building Vcell3D", nthreads())
    @show freeids δ
    filter(net, v) = all(v .!= 0.0) # avoid artifacts at 0.0
    Vcell3D = Dyn.Space(net3D, freeids, δ; filter)
    @sglob Dyn Vcell3D net3D

    ps = Dyn.pol_plots(Vcell3D; bins = 10000)
    sfig(Dyn, ps, "Vcell3D.png")
end


## ------------------------------------------------------
# Test dimensionality
let
    net3D = Dyn.ToyModel3D()
    tol = 1e-3

    Dyn.fix!(net3D, :z, 0.3)
    Dyn.reducebox!(net3D)
    @assert any(net3D.ub .- net3D.lb .> tol)
    
    lb, ub = Dyn.bounds(net3D, :ug)
    Dyn.fix!(net3D, :ug, lb + (ub - lb) / 2)
    Dyn.reducebox!(net3D)
    @assert any(net3D.ub .- net3D.lb .> tol)

    lb, ub = Dyn.bounds(net3D, :uo)
    Dyn.fix!(net3D, :uo, lb + (ub - lb) / 2)
    Dyn.reducebox!(net3D)
    @assert all(net3D.ub .- net3D.lb .≈ 0.0)
    net3D
end