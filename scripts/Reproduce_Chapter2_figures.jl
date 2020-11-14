import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.MonteCarlo
using Plots
using Base.Threads
using ProgressMeter
using Test

## ------------------------------------------------------------------
# Different mrates
@time let
    mutrs = collect(0.0:0.2:1.0)
    ncells, niters = Int(1e7), Int(1e9)
    ps = []
    for mutr in mutrs

        @info "Running Simulation " ncells niters mutr
        
        # simulating 
        p = Polytope()
        pinking_fun(cell) = rand() < vatp(cell)/vatp_global_max(p)
        threading_th = Int(1e5)
        cells = runMC(;p, pinking_fun, ncells, niters, mutr, threading_th, verbose = true)

        # Ploting
        normalize = :probability
        alpha = 0.8
        label = round(mutr, digits = 3)
        plt = histogram(vatp.(cells); normalize, label)

        # Storing
        push!(ps, plt)
        GC.gc()
    end
    plot(ps..., layout = length(ps))
end

## ------------------------------------------------------------------
# Different niters
ps = let
    @info "Different niters"
    niters_order = 5:9
    ncells, mutr = Int(1e7), 0.0
    ps = []
    for order in niters_order
        niters = 10^order
        @info "Running Simulation " ncells niters mutr
        
        # simulating 
        cells = simulation(ncells, niters, mutr, verbose = true)

        # Ploting
        normalize = :probability
        alpha = 0.8
        label = "10^$(order)"
        plt = histogram(vatp.(cells); normalize, label)

        # Storing
        push!(ps, plt)
        GC.gc()
    end
    plot(ps..., layout = length(ps))
end
