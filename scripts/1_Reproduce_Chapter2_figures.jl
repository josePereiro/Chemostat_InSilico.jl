using DrWatson
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.Utilities
using Chemostat_Dynamics.MonteCarlo
using Chemostat_Dynamics.MaxEnt
using Plots
using Statistics
using Base.Threads

## ------------------------------------------------------------------
# TOOLS
function sim_filename(dirtuple; rootdir, name, ext, wkargs...) 
    dir = joinpath(rootdir, savename(dirtuple))
    !isdir(dir) && mkpath(dir)
    joinpath(dir, savename(name, (;dirtuple..., wkargs...), ext))
end

MCdata_filename(dirtuple; wkargs...) = 
    sim_filename(dirtuple; name = "Ch2-MC-Data", rootdir = DATA_DIR, ext = "bson", wkargs...) 
    

MCfig_filename(name, dirtuple; wkargs...) =
    sim_filename(dirtuple; name, rootdir = FIGURES_DIR, ext = "png", wkargs...) 

## ------------------------------------------------------------------
# PLOTS

## ------------------------------------------------------------------
# Figura  2.2:  Comparacion  entre  la  distribucion  de  MaxEnt  (curva  azul)  
# y la  distribucion  estacionaria del modelo estocastico (histograma) para vatpy varios valores de
# Figura 2.4: Distribucion de vglc para diferentes valores demutr.
function plot_f2_2_f2_4(pol, mut_rates, betas, ncells, ciodiff, mstgth)
    
    dirtuple = (;ncells)
    for (fig, marginal_fun, vfun) in [
                                        ("Fig-2.2", vatp_marginal_pdf, vatp), 
                                        ("Fig-2.4", vg_marginal_pdf, vg)
                                    ]
        ps = []
        for (mutr, beta) in zip(mut_rates, betas)

            filename = MCdata_filename(dirtuple; mutr, ciodiff, mstgth)
            MCvatps, MCvgs = load_data(filename)
            cells = Cell.([pol], MCvatps, MCvgs)

            vrange, probs = marginal_fun(pol, beta; n = Int(1e5))   
            
            Δv = step(vrange)
            MaxEnt_v_mean = sum(probs .* vrange .* Δv)
            MC_v_mean = mean(vfun.(cells))
            @info "Ploting $(fig)" mutr beta ciodiff MC_v_mean MaxEnt_v_mean

            normalize = :pdf
            alpha = 0.8
            title = string("mutr: ", round(mutr, digits = 3), " beta: ",  round(beta, digits = 0), 
                "\nciodiff: ", round(ciodiff, digits = 2), " mstgth: ", round(mstgth, digits = 3))
            plt = plot(;title, xlabel = nameof(vfun))
            histogram!(plt, vfun.(cells); normalize, label = "", color = :black)
            plot!(plt, vrange, probs; lw = 3, label = "")
            push!(ps, plt)

            # Saving
            filename = MCfig_filename(fig, dirtuple; mutr, ciodiff, mstgth)
            savefig(plt, filename)
        end
        tplt = plot(ps...; layout = length(ps), titlefont = 9,
                xaxis = nothing, yaxis = nothing, color = :black)
        filename = MCfig_filename("$fig-Total", dirtuple; ciodiff, mstgth)
        savefig(tplt, filename)
    end
end

## ------------------------------------------------------------------
# Figura 2.3: Comparacion entre β y la tasa de mutacion
function plot_f2_3(mut_rates, betas, ncells, ciodiff, mstgth)
    @info "Ploting 2.3"
    title = string("ciodiff: ", round(ciodiff, digits = 2), " mstgth: ", round(mstgth, digits = 3))
    plt = plot(mut_rates, betas; title,
        xlabel = "mut rate", ylabel = "beta", 
        label = "", color = :blue, lw = 3)
    
    dirtuple = (;ncells)
    filename = MCfig_filename("Fig-2.3", dirtuple; ciodiff, mstgth)
    savefig(plt, filename)
end


## ------------------------------------------------------------------
# Figura 2.5: Velocidad media de consumo de glucosa vglc y de secrecion 
# de lactato vlac, inferidas por MaxEnt (ME) comparados con los valores 
# obtenidos de la simulacion en estado estacionario.
function plot_f2_5(pol, mut_rates, betas, ncells, ciodiff, mstgth)
    
    title = string("Correlation mean(vg)", 
        "\nciodiff: ", round(ciodiff, digits = 2), " mstgth: ", round(mstgth, digits = 3))
    vgplt = plot(;title)
    title = string("Correlation mean(vatp)", 
        "\nciodiff: ", round(ciodiff, digits = 2), " mstgth: ", round(mstgth, digits = 3))
    vatpplt = plot(;title)
    dirtuple = (;ncells)

    for (mutr, beta) in zip(mut_rates, betas)

        filename = MCdata_filename(dirtuple; mutr, ciodiff, mstgth)
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

        @info "Ploting 2.5" mutr beta MC_vg_mean MaxEnt_vg_mean MC_vatp_mean MaxEnt_vatp_mean

        scatter!(vgplt, [MaxEnt_vg_mean], [MC_vg_mean]; 
            color = :black, label = "", xlabel = "vatp MaxEnt", ylabel = "vatp MC")
            
        scatter!(vatpplt, [MaxEnt_vatp_mean], [MC_vatp_mean]; 
            color = :black, label = "", xlabel = "vg MaxEnt", ylabel = "vg MC")
    end
    plot!(vgplt, x -> x, 0.09, 0.16; label = "", lw = 3, alpha = 0.5, ls = :dot)
    plot!(vatpplt, x -> x, 1.5, 6.0; label = "", lw = 3, alpha = 0.5, ls = :dot)
    plt = plot([vgplt, vatpplt]..., layout = 2, size = [900, 500])

    dirtuple = (;ncells)
    filename = MCfig_filename("Fig-2.5", dirtuple; ciodiff, mstgth)
    savefig(plt, filename)
end

## ------------------------------------------------------------------
# RUN SIMULATIONS
let
    pol = Polytope()
    mut_rates = 0.0:0.1:1.0 # mutation rate
    mstgths = 0.3:0.1:1.0 # mutation strengths
    nmstgth = 0.0 # not mutation strength
    ncells = Int(5e5)
    @show iters_bkpoints = floor.(Int, ncells .* 10.0 .^ (-1:0.5:3)) |> sort
    curr_bkpoint = 1
    # This must be set so that all breakpoints are possible
    @show feedback_frec = floor(Int, (minimum(iters_bkpoints) ÷ nthreads()) * 0.95) 
    niters = maximum(iters_bkpoints)
    threading_th = Int(1e5)
    gdth = 1e-5
    verbose = true
    dirtuple = (;ncells)
    
    # Monte Carlo
    for mstgth in mstgths, mutr in mut_rates

        # Check files
        all_files_found = true
        for it in iters_bkpoints
            ciodiff = log10(ncells) - log10(it)
            filename = MCdata_filename(dirtuple; mutr, ciodiff, mstgth)
            isfile(filename) ? 
                (@info string(basename(filename), " found")) : 
                (@info string(basename(filename), " missing"); all_files_found = false)
        end
        if all_files_found
            @info "All files founded, skkiping simulation"  mutr
        else
            verbose && @info "Running Montecarlo " ncells mutr niters mstgth

            # Functions
            pinking_fun(cell) = rand() < vatp(cell)/vatpU(pol)
            function feedback_fun(it, cells)
                # Here I take samples at approx iters_bkpoints times
                if curr_bkpoint <= length(iters_bkpoints)
                    bkpoint = iters_bkpoints[curr_bkpoint]
                    if it >= bkpoint 
                        ciodiff = log10(ncells) - log10(bkpoint)
                        filename = MCdata_filename(dirtuple; mutr, ciodiff, mstgth)
                        save_data(filename, (vatp.(cells), vg.(cells)))
                        @assert isfile(filename)
                        curr_bkpoint += 1
                    end
                end
                return false
            end
            
            at_mutate(parent_cell) = generate_similar_cell(parent_cell, rand() * mstgth * Δvatp(pol))
            at_notmutate(parent_cell) = generate_similar_cell(parent_cell, rand() * nmstgth * Δvatp(pol))
            runMC(;p = pol, ncells, niters, mutr, 
                threading_th, pinking_fun, 
                feedback_fun, feedback_frec,
                at_mutate, at_notmutate,
                verbose
            ) 
            curr_bkpoint = 1
        end
    end 

    # Run MaxEnt for each it breakpoint and plot
    for it in iters_bkpoints

        ciodiff = log10(ncells) - log10(it)
        
        # Searching MaxEnt beta
        for mstgth in mstgths
            betas = []
            for mutr in mut_rates
                # Load MC Data
                filename = MCdata_filename(dirtuple; mutr, ciodiff, mstgth)
                MCvatps, MCvgs = load_data(filename)
                cells = Cell.([pol], MCvatps, MCvgs)
                MC_vatp_mean = mean(vatp.(cells))

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
                push!(betas, max(0.0, beta))
            end

            # Ploting
            plot_f2_2_f2_4(pol, mut_rates, betas, ncells, ciodiff, mstgth)
            plot_f2_3(mut_rates, betas, ncells, ciodiff, mstgth)
            plot_f2_5(pol, mut_rates, betas, ncells, ciodiff, mstgth)

        end
    end
    
end