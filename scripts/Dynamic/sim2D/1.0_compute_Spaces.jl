using ProjAssistant
@quickactivate

@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic
end

## ------------------------------------------------------
# common globals
let
    sglob(Dyn;
        δ = 1000
    )
end

## ------------------------------------------------------
@time let
    # input
    net2D = Dyn.ToyModel2D()
    
    freeids = [:z, :ug]
    @lglob Dyn δ
    
    @info("Building Vcell2D")
    @show freeids δ
    filter(net, v) = all(v .!= 0.0) # avoid artifacts at 0.0
    Vcell2D = Dyn.Space(net2D, freeids, δ; filter)
    @sglob Dyn Vcell2D net2D
end