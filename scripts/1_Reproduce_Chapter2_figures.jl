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
data_filename(mutr) = joinpath(DATA_DIR, savename("Ch2-MC-MaxEnt-Results", (;mutr), "bson"))

## ------------------------------------------------------------------
# RUN SIMULATION
pol = Polytope()
mut_rates = 0.0:0.1:1.0
betas = []
for mutr in mut_rates
    # Checking files
    filename = data_filename(mutr)
    if isfile(filename) 
        msg = "Data file already exist, skipping MC. To force recomputation delete the file"
        @info msg filename filesize(filename)
        cells = load_data(filename)
    else
        # Monte Carlo
        ncells, niters = Int(1e6), Int(1e9)
        @info "Running Montecarlo " ncells niters mutr
        pinking_fun(cell) = rand() < vatp(cell)/vatp_global_max(pol)
        threading_th = Int(1e5)
        cells = runMC(;p = pol, pinking_fun, ncells, niters, mutr, 
            threading_th, verbose = true)
        
        save_data(filename, cells)
        @info "Data saved" filename filesize(filename)
    end
    MC_vatp_mean = mean(vatp.(cells))

    # searching MaxEnt beta
    target = [MC_vatp_mean]
    @info "Running Gradient descent " MC_vatp_mean
    x0 = [0.0]
    x1 = [10.0] 
    C = [500.0]
    MaxEnt_vatp_mean = nothing
    beta = grad_desc(;target, x0, x1, C, th = 1e-5, verbose = true) do x
        beta = first(x)
        rvatp, probs = vatp_marginal_pdf(pol, beta; n = Int(1e5))   
        Δvatp = step(rvatp)  
        MaxEnt_vatp_mean = sum(probs .* rvatp .* Δvatp)       
        return MaxEnt_vatp_mean
    end |> first
    push!(betas, beta)

    @info "beta found" beta MC_vatp_mean MaxEnt_vatp_mean

    # Saving
    println(); GC.gc()
end
@assert length(betas) == length(mut_rates)

## ------------------------------------------------------------------
# PLOTS

## ------------------------------------------------------------------
# Figura  2.2:  Comparacion  entre  la  distribucion  de  MaxEnt  (curva  azul)  
# y la  distribucion  estacionaria del modelo estocastico (histograma) para vatpy varios valores de
let
    ps = []
    for (mutr, beta) in zip(mut_rates, betas)
        filename = data_filename(mutr)
        cells = load_data(filename)
        rvatp, probs = vatp_marginal_pdf(pol, beta; n = Int(1e5))   
        
        Δvatp = step(rvatp)
        MaxEnt_vatp_mean = sum(probs .* rvatp .* Δvatp)
        MC_vatp_mean = mean(vatp.(cells))
        @info "Ploting 2.2" mutr beta MC_vatp_mean MaxEnt_vatp_mean

        normalize = :pdf
        alpha = 0.8
        plt = plot(title = string("mutr: ", round(mutr, digits = 3), " beta: ", round(beta, digits = 3)))
        histogram!(plt, vatp.(cells); normalize, label = "", color = :black)
        plot!(plt, rvatp, probs; lw = 3, label = "")
        push!(ps, plt)

        # Saving
        filename = joinpath(FIGURES_DIR, savename("Fig-2.2", (;mutr, beta), "png"))
        savefig(plt, filename)
    end
    tplt = plot(ps...; layout = length(ps), titlefont = 9,
            xaxis = nothing, yaxis = nothing, color = :black)
    filename = joinpath(FIGURES_DIR, "Fig-2.2-Total.png")
    savefig(tplt, filename)
end

## ------------------------------------------------------------------
# Figura 2.3: Comparacion entre β y la tasa de mutacion
let
    @info "Ploting 2.3"
    plt = plot(mut_rates, betas;
        xlabel = "mut rate", ylabel = "beta", 
        label = "", color = :blue, lw = 3)
    filename = joinpath(FIGURES_DIR, "Fig-2.3.png")
    savefig(plt, filename)
end

## ------------------------------------------------------------------
# Figura 2.4: Distribucion de vglc para diferentes valores demutr.
let
    ps = []
    for (mutr, beta) in zip(mut_rates, betas)
        filename = data_filename(mutr)
        cells = load_data(filename)
        rvg, probs = vg_marginal_pdf(pol, beta; n = Int(1e5))   
        
        Δvg = step(rvg)
        MaxEnt_vg_mean = sum(probs .* rvg .* Δvg)
        MC_vg_mean = mean(vg.(cells))
        @info "Ploting 2.4" mutr beta MC_vg_mean MaxEnt_vg_mean

        normalize = :pdf
        alpha = 0.8
        plt = plot(title = string("mutr: ", round(mutr, digits = 3), " beta: ", round(beta, digits = 3)))
        histogram!(plt, vg.(cells); normalize, label = "", color = :black)
        plot!(plt, rvg, probs; lw = 3, label = "")
        push!(ps, plt)

        # Saving
        filename = joinpath(FIGURES_DIR, savename("Fig-2.4", (;mutr, beta), "png"))
        savefig(plt, filename)
    end
    tplt = plot(ps...; layout = length(ps), titlefont = 9,
            xaxis = nothing, yaxis = nothing, color = :black)
    filename = joinpath(FIGURES_DIR, "Fig-2.4-Total.png")
    savefig(tplt, filename)
end

## ------------------------------------------------------------------
# Figura 2.5: Velocidad media de consumo de glucosa vglc y de secrecion 
# de lactato vlac, inferidas por MaxEnt (ME) comparados con los valores 
# obtenidos de la simulacion en estado estacionario.
let
    vgplt, vatpplt = plot(title = "Correlation mean(vg)"), plot(title = "Correlation mean(vatp)")
    for (mutr, beta) in zip(mut_rates, betas)
        filename = data_filename(mutr)
        cells = load_data(filename)
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

    filename = joinpath(FIGURES_DIR, "Fig-2.5.png")
    savefig(plt, filename)
end
