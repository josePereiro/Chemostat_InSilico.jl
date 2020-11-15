using DrWatson
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.Utilities
using Chemostat_Dynamics.MonteCarlo
using Chemostat_Dynamics.MaxEnt
using Plots
using Statistics

## ------------------------------------------------------------------
# TOOLS
MCdata_filename(mutr, ncells, ciodiff) = joinpath(DATA_DIR, 
    savename("Ch2-MC-Data", (;mutr, ncells, ciodiff), "bson"))

## ------------------------------------------------------------------
# SIMULATION FUNCTION
function run_simulation(pol, mut_rates, ncells, ciodiff; 
        threading_th = Int(1e5), gdth = 1e-5, verbose = true)
    betas = []
    niters = floor(Int, ncells .* 10.0^ciodiff)
    for mutr in mut_rates
        filename = MCdata_filename(mutr, ncells, ciodiff)
        if isfile(filename) 
            # Checking files
            msg = "Data file already exist, skipping MC. To force recomputation delete the file"
            verbose && @info msg filename filesize(filename)
            MCvatps, MCvgs = load_data(filename)
            cells = Cell.([pol], MCvatps, MCvgs)
        else
            # Monte Carlo
            verbose && @info "Running Montecarlo " ncells ciodiff mutr niters  
            pinking_fun(cell) = rand() < vatp(cell)/vatp_global_max(pol)
            cells = runMC(;p = pol, pinking_fun, ncells, niters, mutr, 
                    threading_th, verbose)
            
            save_data(filename, (vatp.(cells), vg.(cells)))
            verbose && @info "Data saved" filename filesize(filename)
        end
        MC_vatp_mean = mean(vatp.(cells))

        # searching MaxEnt beta
        target = [MC_vatp_mean]
        verbose && @info "Running Gradient descent " MC_vatp_mean
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
        push!(betas, max(0.0, beta))

        verbose && @info "beta found" beta MC_vatp_mean MaxEnt_vatp_mean
        verbose && println()
        GC.gc()
    end
    @assert length(betas) == length(mut_rates)
    return betas
end

## ------------------------------------------------------------------
# PLOTS

## ------------------------------------------------------------------
# Figura  2.2:  Comparacion  entre  la  distribucion  de  MaxEnt  (curva  azul)  
# y la  distribucion  estacionaria del modelo estocastico (histograma) para vatpy varios valores de
# Figura 2.4: Distribucion de vglc para diferentes valores demutr.
function plot_f2_2_f2_4(pol, mut_rates, betas, ncells, ciodiff)
    for (fig, marginal_fun, vfun) in [
                                        ("Fig-2.2", vatp_marginal_pdf, vatp), 
                                        ("Fig-2.4", vg_marginal_pdf, vg)
                                    ]
        ps = []
        for (mutr, beta) in zip(mut_rates, betas)
            filename = MCdata_filename(mutr, ncells, ciodiff)
            MCvatps, MCvgs = load_data(filename)
            cells = Cell.([pol], MCvatps, MCvgs)

            vrange, probs = marginal_fun(pol, beta; n = Int(1e5))   
            
            Δv = step(vrange)
            MaxEnt_v_mean = sum(probs .* vrange .* Δv)
            MC_v_mean = mean(vfun.(cells))
            @info "Ploting $(fig)" mutr beta ciodiff MC_v_mean MaxEnt_v_mean

            normalize = :pdf
            alpha = 0.8
            title = string("mutr: ", round(mutr, digits = 3), " beta: ", round(beta, digits = 3), "\nciodiff: ", round(ciodiff, digits = 3))
            plt = plot(;title, xlabel = nameof(vfun))
            histogram!(plt, vfun.(cells); normalize, label = "", color = :black)
            plot!(plt, vrange, probs; lw = 3, label = "")
            push!(ps, plt)

            # Saving
            filename = joinpath(FIGURES_DIR, savename(fig, (;mutr, ncells, ciodiff), "png"))
            savefig(plt, filename)
        end
        tplt = plot(ps...; layout = length(ps), titlefont = 9,
                xaxis = nothing, yaxis = nothing, color = :black)
        filename = joinpath(FIGURES_DIR, savename("$(fig)-Total", (;ncells, ciodiff), "png"))
        savefig(tplt, filename)
    end
end

## ------------------------------------------------------------------
# Figura 2.3: Comparacion entre β y la tasa de mutacion
function plot_f2_3(mut_rates, betas, ncells, ciodiff)
    @info "Ploting 2.3"
    title = string("ciodiff: ", round(ciodiff, digits = 3))
    plt = plot(mut_rates, betas; title,
        xlabel = "mut rate", ylabel = "beta", 
        label = "", color = :blue, lw = 3)
        
    filename = joinpath(FIGURES_DIR, savename("Fig-2.3", (;ncells, ciodiff), "png"))
    savefig(plt, filename)
end


## ------------------------------------------------------------------
# Figura 2.5: Velocidad media de consumo de glucosa vglc y de secrecion 
# de lactato vlac, inferidas por MaxEnt (ME) comparados con los valores 
# obtenidos de la simulacion en estado estacionario.
function plot_f2_5(pol, mut_rates, betas, ncells, ciodiff)
    title = string("Correlation mean(vg)\nciodiff: ", round(ciodiff, digits = 3))
    vgplt = plot(;title)
    title = string("Correlation mean(vatp)\nciodiff: ", round(ciodiff, digits = 3))
    vatpplt = plot(;title)
    for (mutr, beta) in zip(mut_rates, betas)
        filename = MCdata_filename(mutr, ncells, ciodiff)
        MCvatps, MCvgs = load_data(filename)
        cells = Cell.([pol], MCvatps, MCvgs)

        rvg, vgprobs = vg_marginal_pdf(pol, beta; n = Int(1e5))   
        Δvg = step(rvg)
        rvatp, vatpprobs = vatp_marginal_pdf(pol, beta; n = Int(1e5))   
        Δvatp = step(rvatp)
        
        MC_vg_mean = mean(vg.(cells))
        MaxEnt_vg_mean = sum(vgprobs .* rvg .* Δvg)
        MC_vatp_mean = mean(vatp.(cells))
        MaxEnt_vatp_mean = sum(vatpprobs .* rvatp .* Δvatp)

        @info "Ploting 1.5" mutr beta MC_vg_mean MaxEnt_vg_mean MC_vatp_mean MaxEnt_vatp_mean

        scatter!(vgplt, [MaxEnt_vg_mean], [MC_vg_mean]; 
            color = :black, label = "", xlabel = "vatp MaxEnt", ylabel = "vatp MC")
            
        scatter!(vatpplt, [MaxEnt_vatp_mean], [MC_vatp_mean]; 
            color = :black, label = "", xlabel = "vg MaxEnt", ylabel = "vg MC")
    end
    plot!(vgplt, x -> x, 0.09, 0.16; label = "", lw = 3, alpha = 0.5, ls = :dot)
    plot!(vatpplt, x -> x, 1.5, 6.0; label = "", lw = 3, alpha = 0.5, ls = :dot)
    plt = plot([vgplt, vatpplt]..., layout = 2, size = [900, 500])

    filename = joinpath(FIGURES_DIR, savename("Fig-2.5", (;ncells, ciodiff), "png"))
    savefig(plt, filename)
end

## ------------------------------------------------------------------
# RUN SIMULATIONS
# Polytope
let
    pol = Polytope()
    mut_rates = 0.0:0.1:1.0
    ncells = Int(1e5)
    ciodiffs = collect(-2:3)
    threading_th = Int(1e5)
    gdth = 1e-5
    for ciodiff in ciodiffs
        betas = run_simulation(pol, mut_rates, ncells, ciodiff; 
            threading_th, gdth, verbose = true)    

        # Ploting
        plot_f2_2_f2_4(pol, mut_rates, betas, ncells, ciodiff)
        plot_f2_3(mut_rates, betas, ncells, ciodiff)
        plot_f2_5(pol, mut_rates, betas, ncells, ciodiff)
    end
end
