using Chemostat_Dynamics
using Chemostat_Dynamics.LP_Implement
using Chemostat_Dynamics.Utilities
using ProgressMeter
using Plots
using Serialization
using Base.Threads
using DrWatson

## ---------------------------------------------------------
function plot_res(M::SimModel; f = (x) -> x)
    p1 = plot(xlabel = "time", ylabel = "conc")
    plot!(p1, f.(M.sg_ts); label = "sg", lw = 3)
    plot!(p1, f.(M.sl_ts); label = "sl", lw = 3)
    
    p2 = plot(xlabel = "time", ylabel = "X")
    plot!(p2, f.(M.X_ts); label = "X", lw = 3)
    
    p = plot([p1, p2]...)
end


## ---------------------------------------------------------
# Find X0
let

    # simModel
    M = SimModel(
        θvatp = 2, 
        θvg = 3, 
        niters = 2000,
        X0 = 1e-5,
        D = 0.001
    )


    run_simulation!(M; verbose = true)

    p = plot_res(M)
    savefig(p, joinpath(CH3_FIGURES_DIR, "test_.png"))
end
##
    # # gradient descent
    # x0 = [1e-2]
    # x1 = [2e-2]
    # C = [1e-3]
    # target = [0.0]
    # X = -1
    # X0 = grad_desc(;target, x0, x1, C, th = 1e-7, verbose = false) do x

    #     local M = SimModel(
    #         θvatp = 2, 
    #         θvg = 3, 
    #         niters = 10,
    #         X0 = first(x)
    #     )

    #     lX = last(run_simulation!(M; cache, verbose = true).X_ts)
    #     dX = abs(lX - X) / max(lX, X)
    #     X = lX
    #     return dX
    # end |> first

# end

## ---------------------------------------------------------
target = [MC_vatp_mean]
                verbose && @info "Running Gradient descent " mutr ciodiff mstgth MC_vatp_mean
                x0 = [0.0]
                x1 = [10.0] 
                C = [500.0]
                MaxEnt_vatp_mean = nothing
                beta = grad_desc(;target, x0, x1, C, th = gdth, verbose) do x
                    beta = first(x)
                    rvatp, probs = vatp_marginal_pdf(pol, beta; n = Int(1e5))   
                    Δvatp = step(rvatp)  
                    MaxEnt_vatp_mean = sum(probs .* rvatp .* Δvatp)       
                    return MaxEnt_vatp_mean
                end |> first
## ---------------------------------------------------------
# # ToyModel Simulations
# let
#     lagtime = 1000
#     etime = 5000
#     Ds = collect(10.0.^(-7:0.5:-2))
#     bkps = collect(lagtime .+ etime * (1:length(Ds)))
#     bkpi = 1

#     # cache
    

#     P = SimParams(;
#         θvatp = 2, θvg = 3,
#         D = Ds[bkpi], 
#         X0 = 1e-6, 
#         Xmin = 1e-20, 
#         niters = lagtime + etime * length(Ds) + etime,
#         damp = 0.95,
#         ϵ = 0.0,
#         sg0 = 5.0,
#         cache_marginf = 0.1
#     )

#     function at_iter(it, P) 
#         next_bkp = bkps[bkpi]
#         (it != next_bkp || bkpi == length(Ds)) && return P
#         bkpi += 1
#         new_D = Ds[bkpi]
#         return SimParams(P; D = new_D)
#     end
#     X_ts, sg_ts, sl_ts, Xb = run_simulation(P; at_iter)

#     fname = mysavename("Ch3sim_"; P.X0, P.niters, P.ϵ)
#     @show fname
#     simfile = joinpath(CH3_DATA_DIR, string(fname, ".jld"))
#     serialize(simfile, (X_ts, sg_ts, sl_ts, Xb, P, Ds))
#     # X_ts, sg_ts, sl_ts, Xb, P, Ds = deserialize(simfile)
#     p1 = plot(xlabel = "time", ylabel = "conc")
#     plot!(p1, log10.(sg_ts); label = "sg", lw = 3)
#     plot!(p1, log10.(sl_ts); label = "sl", lw = 3)
    
#     p2 = plot(xlabel = "time", ylabel = "X")
#     plot!(p2, log10.(X_ts); label = "X", lw = 3)
    
#     p = plot([p1, p2]...)

#     # save fig
#     figname = joinpath(CH3_FIGURES_DIR, string(fname, ".png"))
#     savefig(p, figname)
# end
