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
@time let
    # input
    net = Dyn.ToyModel2D()
    Dyn.reducebox!(net) # This is important
    
    freeids = [:z, :ug]
    freeis = Dyn.rxnindex(net, freeids)
    stoitolf = 0.01
    δ = 1000
    
    # container, net, left_freeis, head_freeis, ϵ
    @info("Building Vcell0")
    @show freeids stoitolf δ
    Vcell0 = Dyn.Space(net, freeids, δ; stoitolf)
    Dyn.sglob(;net, Vcell0, stoitolf, freeids, freeis, δ)

end;

## ------------------------------------------------------
let
    Vcell0 = Dyn.lglob(:Vcell0)
    p = Dyn.plot_marginals(Vcell0; normalize = false)
end