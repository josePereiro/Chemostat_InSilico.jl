import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.MonteCarlo
using Plots

## ------------------------------------------------------------------
# Plot, Sampling histogram
@time let
    n = 1_000_000
    p = Polytope()
    cells_pool = generate_random_cells(p, n; verbose = false)
    mvatp = vatpU(p)
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
@time let
    n = 10_000
    p = Polytope()
    cells_pool = generate_random_cells(p, n; verbose = false)
    mvatp = vatpU(p)
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
@time let
    pol = Polytope(xi = 100.0)
    betas_orders = -5:0.5:5
    plt = plot()
    for o in betas_orders
        rvatp, probs = vatp_marginal_pdf_probs(pol, 10.0^o; n = Int(1e5))
        plot!(plt, rvatp, probs ./ maximum(probs), label = "", lw = 3)
    end
    plt
end

## ------------------------------------------------------------------
# Different mrates
@time let
    mutrs = collect(0.0:0.2:1.0)
    ncells, niters = Int(1e6), Int(1e7)
    ps = []
    for mutr in mutrs

        @info "Running Simulation " ncells niters mutr
        
        # simulating 
        p = Polytope()
        pinking_fun(cell) = rand() < vatp(cell)/vatpU(p)
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
@time let
    @info "Different niters"
    niters_order = 5:8
    ncells, mutr = Int(1e6), 0.0
    pol = Polytope()
    ps = []
    for order in niters_order
        niters = 10^order
        @info "Running Simulation " ncells niters mutr
        
        # simulating 
        pinking_fun(cell) = rand() < vatp(cell)/vatpU(pol)
        threading_th = Int(1e5)
        cells = runMC(;p = pol, pinking_fun, ncells, niters, mutr, threading_th, verbose = true)

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
