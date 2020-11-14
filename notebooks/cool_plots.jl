import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import Chemostat_Dynamics: SimpleDynamic
const SD = SimpleDynamic
using Plots

## ------------------------------------------------------------------
# Plot, Sampling histogram
let
    n = 1_000_000
    p = SD.SimParams()
    cells_pool = SD.generate_random_cells(p, n; verbose = false)
    mvatp = SD.vatp_global_max(p)
    pcells = SD.pick_cells(n, cells_pool) do cell
        prob = SD.vatp(cell)/mvatp
        return rand() <= prob
    end
    
    # Ploting
    ps = []
    for (name, f) in [("vatp", SD.vatp), ("vg", SD.vg)]
        plt = plot(title = "Polytope Sampling",legend = :left, xlabel = name, ylabel = "prob")
        normalize = :probability
        alpha = 0.8
        Plots.histogram!(plt, f.(cells_pool); label = "uniform", normalize, alpha)
        Plots.histogram!(plt, f.(pcells); label = "non-uniform", normalize, alpha)
        push!(ps, plt)
    end
    plot(ps...; layout = 2, size = [850, 400])
end

## ------------------------------------------------------------------
# Plot, Sampling polytope graphs
let
    n = 10_000
    p = SD.SimParams()
    cells_pool = SD.generate_random_cells(p, n; tries = 100, verbose = true)
    mvatp = SD.vatp_global_max(p)
    pcells = SD.pick_cells(n, cells_pool) do cell
        prob = SD.vatp(cell)/mvatp
        return rand() <= prob
    end
    ps = []
    for (title, cells) in [("Uniform sampling", cells_pool), 
                           ("Non-Uniform Sampling", pcells)]                
        plt = SD.plot_polytop(p)
        plot!(plt; title)
        scatter!(plt, SD.vatp.(cells), SD.vg.(cells), label = "", color = :black, alpha = 0.1)
        push!(ps, plt)
    end
    plot(ps...; layout = 2, size = [850, 400])
end


