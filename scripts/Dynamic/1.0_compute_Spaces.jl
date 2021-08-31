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
        # δ = 1000
        δ = 200
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
 
## ------------------------------------------------------
@time let
    # input
    net3D = Dyn.ToyModel3D()
    
    freeids = [:z, :ug, :uo]
    @lglob Dyn δ
    
    @info("Building Vcell3D")
    @show freeids δ
    filter(net, v) = all(v .!= 0.0) # avoid artifacts at 0.0
    Vcell3D = Dyn.Space(net3D, freeids, δ; filter)
    @sglob Dyn Vcell3D net3D
end

# ## ------------------------------------------------------
# # dev
# using Plots
# let
#     net = Dyn.ToyModel3D()
 
#     # Dyn.fix!(net, :z, 8.0)
#     # Dyn.fix!(net, :ug, 15.0)
#     # Dyn.fix!(net, :uo, 4.0)
#     # Dyn.reducebox!(net)
#     L, U = Dyn.fva(net)
    
#     @show net
#     @show L
#     @show U
    
#     p = plot()
#     scatter!(p, net.rxns, L; label = "", c = :blue)
#     plot!(p, net.rxns, L; label = "L", c = :blue, ls = :dash)
#     scatter!(p, net.rxns, U; label = "", c = :red)
#     plot!(p, net.rxns, U; label = "U", c = :red, ls = :dash)
 
#     sfig(Dyn, p, "test.png")
#  end
