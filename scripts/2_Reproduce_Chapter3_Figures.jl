import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

import Chemostat_InSilico
const InLP = Chemostat_InSilico.LP_Implement
const InU = Chemostat_InSilico.Utilities

using ProgressMeter
using Plots
using Serialization
using Base.Threads
using Dates
using Base.Threads
using BenchmarkTools

## ---------------------------------------------------------
# let
#     M = InLP.SimModel(;
#             θvatp = 2, 
#             θvg = 3
#         )

#     # InLP.plot_politope(M)
#     InLP.plot_res(M, InLP.ResTS())
# end

## ---------------------------------------------------------
# Find X0
let

    sim_name = "Ch3_$(now())"
    @info "Doing" sim_name

    M0 = InLP.SimModel(;
            θvatp = 2, 
            θvg = 3, 
            niters = Int(1e8),
            sg0 = 4.5,
            X0 = 0.3,
            damp = 0.0,
            D = 0.01
        )
        
    # cache
    InLP.vgvatp_cache(M0)

    # simulation
    M = deepcopy(M0)

    
    # TODO: make a real chemostat. move D to get a given X
    save_frec = 100
    write_frec = 10
    ts = InLP.ResTS()
    function at_iter(it, M)
        if (it % write_frec == 0)
            push!(ts, M)
        end

        # save fig
        if (it == 1 || it % save_frec == 0 || it == M.niters)

            fname = InLP.mysavename(sim_name, "png"; it, M.D, M.damp)
            path = joinpath(InLP.CH3_FIGURES_DIR, fname)
            p = InLP.plot_res(M, ts)
            savefig(p, path)

        end
    end
    InLP.run_simulation!(M; at_iter, verbose = true)

end
