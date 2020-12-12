import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

using Chemostat_InSilico
using Chemostat_InSilico.LP_Implement
using Chemostat_InSilico.Utilities
using ProgressMeter
using Plots
using Serialization
using Base.Threads

## ---------------------------------------------------------
using BenchmarkTools
let
    N = 5000
    intr = Int.(1:N)
    strr = string.(1:N)
    d1 = Dict{Tuple{Int64,String}, String}((i, s) => "bla" for i in intr, s in strr)
    d2 = Dict{Int64, Dict{String, String}}(i => Dict{String, String}(s => "bla" for s in strr) for i in intr)
    @assert length(d1) == sum(length.(values(d2)))

    @info "Accessing d1"
    @btime begin 
        for i in $intr, s in $strr
            $d1[(i, s)] = "blo"
        end
    end

    @info "Accessing d2 no-caching"
    @btime begin 
        for i in $intr, s in $strr
            $d2[i][s] = "blo"
        end
    end

    @info "Accessing d2 caching"
    @btime begin 
        for i in $intr
            ld = $d2[i]
            for s in $strr
                ld[s] = "blo"
            end
        end
    end
    @assert length(d1) == sum(length.(values(d2)))
end

##
    # @btime begin 
    #     for _ in 1:100
    #         for i in 1:N
    #             ld::Dict{Int, String} = $d2[i] 
    #             for s in 1:N
    #                 ld[Int(s)] = "blo"
    #             end
    #         end
    #     end
    # end

    

# end

## ---------------------------------------------------------
# Find X0
let

    M0 = SimModel(;
            θvatp = 2, 
            θvg = 3, 
            niters = 1000000,
            sg0 = 4.5,
            X0 = 0.3,
            damp = 0.98
        )
        
    writing_lock = ReentrantLock()
        
    # cache
    vgvatp_cache(M0)
        
    # Ds = collect(10.0.^(-3:0.1:-1))
    # @threads for D in Ds

    #     # simModel
    #     M = deepcopy(M0)
    #     M.D = D

    #     save_frec = M.niters ÷ 100

    #     lock(writing_lock) do
    #         @info "Doing" threadid() D
    #         println()
    #     end

    #     # TODO: make a real chemostat. move D to get a given X
    #     function at_iter(it, M)
    #         (it == 1 || it % save_frec != 0 || it == M.niters) && return

    #         lock(writing_lock) do
    #             @info "Report" threadid() it
    #             println()
    #         end    
    #     end
    #     run_simulation!(M; verbose = false)

    #     datname = mysavename("Ch3_", "jld"; M.D, M.damp)
    #     serialize(joinpath(CH3_DATA_DIR, datname), M)

    #     lock(writing_lock) do
    #         @info "Finised" threadid() D
    #         println()
    #     end    
    # end

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
# target = [MC_vatp_mean]
#                 verbose && @info "Running Gradient descent " mutr ciodiff mstgth MC_vatp_mean
#                 x0 = [0.0]
#                 x1 = [10.0] 
#                 C = [500.0]
#                 MaxEnt_vatp_mean = nothing
#                 beta = grad_desc(;target, x0, x1, C, th = gdth, verbose) do x
#                     beta = first(x)
#                     rvatp, probs = vatp_marginal_pdf(pol, beta; n = Int(1e5))   
#                     Δvatp = step(rvatp)  
#                     MaxEnt_vatp_mean = sum(probs .* rvatp .* Δvatp)       
#                     return MaxEnt_vatp_mean
#                 end |> first
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
