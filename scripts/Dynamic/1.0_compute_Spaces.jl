@time begin 
    import Chemostat_InSilico
    const Dyn = Chemostat_InSilico.Dynamic

    import UtilsJL
    const PltU = UtilsJL.PlotsUtils
    const Ass = UtilsJL.ProjAssistant
    const GU = UtilsJL.GeneralUtils

    Ass.set_verbose(false)
end

## ------------------------------------------------------
# common globals
let
    Dyn.sglob(;
        δ = 1000
    )
end

## ------------------------------------------------------
@time let
    # input
    net2D = Dyn.ToyModel2D()
    
    freeids = [:z, :ug]
    δ = Dyn.lglob(:δ)
    
    @info("Building Vcell2D")
    @show freeids δ
    filter(net, v) = all(v .!= 0.0) # avoid artifacts at 0.0
    Vcell2D = Dyn.Space(net2D, freeids, δ; filter)
    Dyn.sglob(;Vcell2D, net2D)

end

## ------------------------------------------------------
@time let
    # input
    net3D = Dyn.ToyModel3D()
    
    freeids = [:z, :ug, :uo]
    # δ = Dyn.lglob(:δ)
    δ = 1000
    
    @info("Building Vcell3D")
    @show freeids δ
    filter(net, v) = all(v .!= 0.0) # avoid artifacts at 0.0
    Vcell3D = Dyn.Space(net3D, freeids, δ; filter)
    Dyn.sglob(;Vcell3D, net3D)
end