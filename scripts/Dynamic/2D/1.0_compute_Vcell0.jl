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
        stoitolf = 0.01, 
        δ = 1000
    )
end

## ------------------------------------------------------
@time let
    # input
    net = Dyn.ToyModel2D()
    Dyn.reducebox!(net) # This is important
    
    freeids = [:z, :ug]
    freeis = Dyn.rxnindex(net, freeids)
    stoitolf = Dyn.lglob(:stoitolf)
    δ = Dyn.lglob(:δ)
    
    # container, net, left_freeis, head_freeis, ϵ
    @info("Building Vcell2D")
    @show freeids stoitolf δ
    Vcell2D = Dyn.Space(net, freeids, δ; stoitolf)
    Dyn.sglob(;Vcell2D)

    Dyn.plot_marginals(Vcell2D; normalize = false)
end;

## ------------------------------------------------------
@time let
    # input
    net = Dyn.ToyModel3D()
    Dyn.reducebox!(net) # This is important
    
    freeids = [:z, :ug, :uo]
    freeis = Dyn.rxnindex(net, freeids)
    stoitolf = Dyn.lglob(:stoitolf)
    δ = Dyn.lglob(:δ)
    
    # container, net, left_freeis, head_freeis, ϵ
    @info("Building Vcell3D")
    @show freeids stoitolf δ
    Vcell3D = Dyn.Space(net, freeids, δ; stoitolf)
    Dyn.sglob(;Vcell3D)
end;
