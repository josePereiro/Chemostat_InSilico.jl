import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.MonteCarlo
using Plots

## ------------------------------------------------------------------
# Plot, Sampling histogram
let
    n = 1_000_000
    p = Polytope()
    cells_pool = generate_random_cells(p, n; verbose = false)
    mvatp = vatp_global_max(p)
    pcells = pick_cells(n, cells_pool) do cell
        prob = vatp(cell)/mvatp
        return rand() <= prob
    end
    
    # Ploting
    ps = []
    for (name, f) in [("vatp", vatp), ("vg", vg)]
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
    p = Polytope()
    cells_pool = generate_random_cells(p, n; verbose = false)
    mvatp = vatp_global_max(p)
    pcells = pick_cells(n, cells_pool) do cell
        prob = vatp(cell)/mvatp
        return rand() <= prob
    end
    ps = []
    for (title, cells) in [("Uniform sampling", cells_pool), 
                           ("Non-Uniform Sampling", pcells)]                
        plt = plot_polytope(p)
        plot!(plt; title)
        scatter!(plt, vatp.(cells), vg.(cells), label = "", color = :black, alpha = 0.1)
        push!(ps, plt)
    end
    plot(ps...; layout = 2, size = [850, 400])
end


## ------------------------------------------------------------------
# Plot, vatp vs beta
let
    pol = Polytope(xi = 100.0)
    betas_orders = -5:0.5:5
    plt = plot()
    for o in betas_orders
        rvatp, probs = vatp_marginal_probs(pol, 10.0^o; n = Int(1e5))
        plot!(plt, rvatp, probs ./ maximum(probs), label = "", lw = 3)
    end
    plt
end
