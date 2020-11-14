import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import Chemostat_Dynamics: SimpleDynamic
const SD = SimpleDynamic
using Base.Threads
using ProgressMeter
using Plots
using Test

## ------------------------------------------------------------------
function generate_dougther(c)
    vglb, vgub = SD.vg_local_min(c.vatp, c.p), SD.vg_local_max(c.vatp, c.p)
    rvg = rand()*(vgub - vglb) + vglb
    SD.Cell(c.p, c.vatp, rvg)
end

## ------------------------------------------------------------------
function simulation(ncells, niters, mutr; verbose = true)
    # Params and tools
    p = SD.SimParams()
    mvatp = SD.vatp_global_max(p)
    pfun(cell) = rand() < SD.vatp(cell)/mvatp
    itmutates() = rand() < mutr

    # Generating pool
    cells_pool = SD.generate_random_cells(p, ncells; verbose)
    
    # Simulating time
    verbose && (prog = Progress(niters, "Simulating... "))
    for it in 1:niters
        pcell = first(SD.pick_cells(pfun, 1, cells_pool))
        ridx = rand(1:ncells)
        cells_pool[ridx] = itmutates() ? first(SD.generate_random_cells(p, 1)) : generate_dougther(pcell)
        verbose && next!(prog)
    end
    verbose && finish!(prog)
    cells_pool
end

## ------------------------------------------------------------------
ps = let
    mutrs = collect(0.0:0.2:1.0)
    ncells, niters = Int(1e6), Int(1e6)
    ps = []
    for mutr in mutrs

        @info "Running Simulation " ncells niters mutr
        
        # simulating 
        cells = simulation(ncells, niters, mutr, verbose = true)

        # Ploting
        normalize = :probability
        alpha = 0.8
        label = round(mutr, digits = 3)
        plt = histogram(SD.vatp.(cells); normalize, label)

        # Storing
        push!(ps, plt)
        GC.gc()
    end
    plot(ps..., layout = length(ps))
end
