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

# ---------------------------------------------------------
# let
#     N = 100
#     chuncks = InU.get_chuncks(1:N, 4, th = -1)
#     for _ in 1:1000
#         d = Dict()
#         for chunck in chuncks, i in chunck
#             d[i] = rand()
#         end
#         @threads for chunck in chuncks
#             for i in chunck
#                 d[i] = rand()
#             end
#         end
#         length(d) |> println
#     end
# end
# ---------------------------------------------------------
let
    M0 = InLP.SimModel(;
        θvatp = 2, 
        θvg = 3, 
        niters = Int(1e8),
        sg0 = 4.5,
        X0 = 0.3,
        damp = 0.98,
        D = 0.01
    )
    vatp_range, vg_ranges = InLP.vatpvg_ranges(M0)
    InU.get_chuncks(vatp_range, 4, th = 1)
end
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
            damp = 0.98,
            D = 0.01
        )
        
    # cache
    InLP.vgvatp_cache(M0)

    # simulation
    M = deepcopy(M0)

    save_frec = 1000

    # TODO: make a real chemostat. move D to get a given X
    function at_iter(it, M)
        if (it == 1 || it % save_frec == 0 || it == M.niters)

            fname = InLP.mysavename(sim_name, "png"; it, M.D, M.damp)
            # serialize(joinpath(InLP.CH3_DATA_DIR, datname), M)
            path = joinpath(InLP.CH3_FIGURES_DIR, fname)
            p = InLP.plot_res(M)
            savefig(p, path)

        end
    end
    InLP.run_simulation!(M; at_iter, verbose = true)

end
##
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
    #     InLP.run_simulation!(M; verbose = false)

    #     datname = InLP.mysavename("Ch3_", "jld"; M.D, M.damp)
    #     serialize(joinpath(InLP.CH3_DATA_DIR, datname), M)

    #     lock(writing_lock) do
    #         @info "Finised" threadid() D
    #         println()
    #     end    
    # end

# end

##
    # # gradient descent
    # x0 = [1e-2]
    # x1 = [2e-2]
    # C = [1e-3]
    # target = [0.0]
    # X = -1
    # X0 = grad_desc(;target, x0, x1, C, th = 1e-7, verbose = false) do x

    #     local M = InLP.SimModel(
    #         θvatp = 2, 
    #         θvg = 3, 
    #         niters = 10,
    #         X0 = first(x)
    #     )

    #     lX = last(InLP.run_simulation!(M; cache, verbose = true).X_ts)
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

#     fname = InLP.mysavename("Ch3sim_"; P.X0, P.niters, P.ϵ)
#     @show fname
#     simfile = joinpath(InLP.CH3_DATA_DIR, string(fname, ".jld"))
#     serialize(simfile, (X_ts, sg_ts, sl_ts, Xb, P, Ds))
#     # X_ts, sg_ts, sl_ts, Xb, P, Ds = deserialize(simfile)
#     p1 = plot(xlabel = "time", ylabel = "conc")
#     plot!(p1, log10.(sg_ts); label = "sg", lw = 3)
#     plot!(p1, log10.(sl_ts); label = "sl", lw = 3)
    
#     p2 = plot(xlabel = "time", ylabel = "X")
#     plot!(p2, log10.(X_ts); label = "X", lw = 3)
    
#     p = plot([p1, p2]...)

#     # save fig
#     figname = joinpath(InLP.CH3_FIGURES_DIR, string(fname, ".png"))
#     savefig(p, figname)
# end

# ## ---------------------------------------------------------
# using BenchmarkTools
# let
#     N = 5000
#     intr = Int.(1:N)
#     strr = string.(1:N)
#     d1 = Dict{Tuple{Int64,String}, String}((i, s) => "bla" for i in intr, s in strr)
#     d2 = Dict{Int64, Dict{String, String}}(i => Dict{String, String}(s => "bla" for s in strr) for i in intr)
#     @assert length(d1) == sum(length.(values(d2)))

#     @info "Accessing d1"
#     @btime begin 
#         for i in $intr, s in $strr
#             $d1[(i, s)] = "blo"
#         end
#     end

#     @info "Accessing d2 no-caching"
#     @btime begin 
#         for i in $intr, s in $strr
#             $d2[i][s] = "blo"
#         end
#     end

#     @info "Accessing d2 caching"
#     @btime begin 
#         for i in $intr
#             ld = $d2[i]
#             for s in $strr
#                 ld[s] = "blo"
#             end
#         end
#     end
#     @assert length(d1) == sum(length.(values(d2)))
# end

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
# using BenchmarkTools
# let
#     N, M = 1000, 1600
#     # is = 1:rand(1:N)
#     is = 1:N
#     # js = [1:rand(1:M) for i in is]
#     js = [1:M for i in is]

#     function fun1(d, is, js)
#         unfea_is = setdiff(keys(d), is)
#         foreach((i) -> delete!(d, i), unfea_is)
#         for i in is
#             ld = d[i]
#             unfea_js = setdiff(keys(ld), js[i])
#             foreach((j) -> delete!(ld, j), unfea_js)
#         end
#     end

#     function fun2(d, is, js)
#         unfea_is = setdiff(keys(d), is)
#         for i in is
#             ld = d[i]
#             unfea_js = setdiff(keys(ld), js[i])
#         end
#     end

#     function fun3(d, is, js)
#         for (i, di) in d
#             foreach((j) -> di[j] = 0.0, keys(di))
#         end

#         for i in is
#             di = d[i]
#             @inbounds for j in js[i]
#                 di[j] = 1.0
#             end
#         end
#     end

#     # @info "fun1"
#     # d1 = Dict(i => Dict(j => rand() for j in 1:M) for i in 1:N)
#     # @btime $fun1($d1, $is, $js)

#     # @info "fun2"
#     # d2 = Dict(i => Dict(j => rand() for j in 1:M) for i in 1:N)
#     # @btime $fun2($d2, $is, $js)

#     @info "fun3"
#     d3 = Dict(i => Dict(j => rand() for j in 1:M) for i in 1:N)
#     @btime $fun3($d3, $is, $js)

#     # @assert sum(length.(values(d1))) == sum(length.(values(d3)))
# end

# ## ---------------------------------------------------------
# let
#     M, N = 3000, 3000
#     d = Dict{Int, Dict{Int, Float64}}(i => Dict(j => rand() for j in 1:M) for i in 1:N)

#     function fun1(d::Dict{Int, Dict{Int, Float64}})
#         c = 0.0
#         for i in 1:N
#             di = get!(d, i, Dict{Int, Float64}())
#             for j in 1:M
#                 c += get!(di, j, 0.0)
#             end
#         end
#         c
#     end

#     function fun2(d::Dict{Int, Dict{Int, Float64}})
#         c = 0.0
#         for i in 1:N
#             di = haskey(d, i) ? d[i] : Dict{Int, Float64}()
#             for j in 1:M
#                 c += haskey(di, j) ? di[j] : 0.0
#             end
#         end
#         c
#     end

#     function fun3(d::Dict{A, Dict{B, C}}) where {A, B, C}
#         c = zero(C)
#         for i in 1:N
#             di = haskey(d, i) ? d[i] : Dict{B, C}()
#             for j in 1:M
#                 c += haskey(di, j) ? di[j] : zero(C)
#             end
#         end
#         c
#     end

#     # function fun3(d::Dict{Int, Dict{Int, Float64}})
#     #     c = 0.0
#     #     for i in 1:N
#     #         di::Dict{Int, Float64} = if haskey(d, i); d[i]; else Dict{Int, Float64}() end
#     #         for j in 1:M
#     #             c += if haskey(di, j); di[j]; else 0.0 end::Float64
#     #         end
#     #     end
#     #     c
#     # end

#     d1 = deepcopy(d)
#     @info "fun1" fun1(d1)
#     @btime $fun1($d1)
    
#     d2 = deepcopy(d)
#     @info "fun2" fun2(d2)
#     @btime $fun2($d2)
    
#     d3 = deepcopy(d)
#     @info "fun3" fun3(d3)
#     @btime $fun3($d3)

#     # @code_warntype fun2(d)

#     nothing
# end