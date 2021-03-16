import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin
    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    using Base.Threads
    using Dates
    using Statistics
    using InteractiveUtils
    using Random

    using Plots
    using ProgressMeter
    using Base.Threads
    using ColorSchemes
    UJL.set_cache_dir(Dyn.cachedir())
    UJL.set_verbose(false)
end

## -----------------------------------------------------------------------------------------------
fileid = "5"
mysavefig(p, pname; params...) = 
    Dyn.mysavefig(p, pname, Dyn.plotsdir(), fileid; params...)

function entropy(P)
    H = 0.0
    for (vatp, mP) in P
        for (vg, p) in mP
            iszero(p) && continue
            H -= p * log2(p)
        end
    end
    @assert !isnan(H)
    H
end

getbeta(order) = sign(order) * 10.0 ^ abs(order)

function get_ticks(f::Function, values; l = length(values))
    idxs = unique(floor.(Int, range(1, length(values); length = l)))
    values = values[idxs]
    values, f.(values)
end

## -----------------------------------------------------------------------------------------------
# Compute beta maps
BEXP_DID = (:BETA_SPACE_EXPLORATION_RESULTS)

DTASK = @async let
    error("Comment this line to continue")
    # base model
    M = Dyn.SimModel(;
        δvatp = 2, 
        δvg = 3, 
    )

    biom_idx = M.obj_idx
    LP_cache = Dyn.vgvatp_cache(M)
    z(vatp, vg) = LP_cache[vatp][vg][biom_idx]
    vg(vatp, vg) = vg

    beta_orders = -3:0.1:3

    H_board = zeros(length(beta_orders), length(beta_orders))
    vg_board = similar(H_board)
    z_board = similar(H_board)
    prog = Progress(length(H_board))

    Ch = Channel(nthreads()) do Ch_
        for (b1i, order1) in enumerate(beta_orders)
            beta1 = getbeta(order1)
            for (b2i, order2) in enumerate(beta_orders)
                beta2 = getbeta(order2)
                put!(Ch_, (b1i, beta1, b2i, beta2))
            end
        end
    end
    
    @threads for _ in 1:nthreads()
        thid = threadid()
        thid == 1 && continue # Let main free
        
        local M0 = deepcopy(M)
        local PME = Dyn.get_join(M0)
        for (b1i, beta1, b2i, beta2) in Ch
            
            try;
                PME = Dyn.get_join!(M0, PME) do _vatp, _vg
                    exp(beta1 * z(_vatp, _vg) + beta2 * vg(_vatp, _vg))
                end
                H = entropy(PME)

                H_board[b1i, b2i] = H
                vg_board[b1i, b2i] = Dyn.ave_over(PME) do vatp, vg
                    vg
                end
                z_board[b1i, b2i] = Dyn.ave_over(PME) do vatp, vg
                    z(vatp, vg)
                end

                # next!(prog; showvalues = [
                #         (:z_beta, beta1),    
                #         (:vg_beta, beta2),    
                #         (:H, H),    
                #         (:thid, thid),    
                #     ]
                # )
            catch err
                @warn("Error at", beta1, beta2)
                rethrow(err)
            end
        end

    end # for thid

    # save
    dat = (;M, beta_orders, H_board, vg_board, z_board)
    UJL.save_cache(BEXP_DID, dat; verbose = true)

end 

## -----------------------------------------------------------------------------------------------
# raw plots
let
    M, beta_orders, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)

    ticks = get_ticks(beta_orders; l = 11) do order
        round(getbeta(order))
    end

    common_params = (;
        xticks = ticks,
        yticks = ticks,
        xrotation = 35,
    )
    p1 = heatmap(beta_orders, beta_orders, H_board; 
        title = "Entropy", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    )
    p2 = heatmap(beta_orders, beta_orders, vg_board; 
        title = "Mean vg", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    )
    p3 = heatmap(beta_orders, beta_orders, z_board; 
        title = "Mean z", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    )
    ps = Plots.Plot[p1, p2, p3]
    
    # saving
    layout = (1, 3)
    mysavefig(ps, "raw_maps")
end

## -----------------------------------------------------------------------------------------------
# entropy plots
let
    M, beta_orders, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    betas = getbeta.(beta_orders)
    maxbeta = maximum(betas)

    xticks = get_ticks(eachindex(beta_orders); l = 11) do idx
        round(getbeta(beta_orders[idx]))
    end
    xrotation = 35
    colors = palette([:blue, :red], length(beta_orders))

    # H_board[zb, vgb]
    ps = Plots.Plot[]
    lw = 3
    alpha = 0.8
    p = plot(;xlabel = "vg beta", ylabel = "H", xticks, xrotation)
    for (color, row) in zip(colors, eachrow(H_board))
        p = plot!(p, row; label = "", lw, alpha, color)
    end
    push!(ps, p)

    p = plot(;xlabel = "z beta", ylabel = "H", xticks, xrotation)
    for (color, col) in zip(colors, eachcol(H_board))
        p = plot!(p, col; label = "", lw, alpha, color)
    end
    push!(ps, p)

    mysavefig(ps, "entropy_vs_beta")
end

## -----------------------------------------------------------------------------------------------
# beta mag vs entropy
let
    M, beta_orders, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    betas = getbeta.(beta_orders)

    mag(vs) = sqrt(sum(vs .^ 2))

    mags = []
    Hs = []
    for (b1i, beta1) in betas |> enumerate
        for (b2i, beta2) in betas |> enumerate
            push!(mags, mag([beta1, beta2]))
            push!(Hs, H_board[b1i, b2i])
        end
    end

    p = scatter(mags, Hs; 
        label = "", xlabel = "beta mag", ylabel = "entropy", color = :black
    )
    mysavefig(p, "beta_mag_vs_entropy")
end

## -----------------------------------------------------------------------------------------------
function get_validity_boards(M, betas, z_board, vg_board; z_th = 0.05)
    cgD_X = M.cg * M.D / M.X
    z_valid_board = similar(z_board, Bool)
    vg_valid_board = similar(z_board, Bool)
    z_th = M.D * z_th
    for (b1i, beta1) in betas |> enumerate
        for (b2i, beta2) in betas |> enumerate
            z = z_board[b1i, b2i]
            vg = vg_board[b1i, b2i]

            # check validity
            z_valid_board[b1i, b2i] = abs(z - M.D) < z_th
            vg_valid_board[b1i, b2i] = vg < cgD_X
        end
    end
    z_valid_board, vg_valid_board
end

## -----------------------------------------------------------------------------------------------
# validity plots
let
    M, beta_orders, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    betas = getbeta.(beta_orders)

    # fake stst params
    cgD_X = 0.3
    M.D = 0.03
    M.X = 1.0
    M.cg = M.X * cgD_X / M.D
    
    z_th = 0.08
    z_valid_board, vg_valid_board = get_validity_boards(M, betas, z_board, vg_board; z_th)

    ticks = get_ticks(beta_orders; l = 11) do order
        round(getbeta(order))
    end

    common_params = (;
        xticks = ticks,
        yticks = ticks,
        xrotation = 35,
        grid = false,
    )
    
    ps = Plots.Plot[]

    # validity maps
    function just_valids(dat, validity)
        vals = dat .* validity
        ifelse.(iszero.(vals), NaN, vals)
    end

    p = heatmap(beta_orders, beta_orders, just_valids(z_board, z_valid_board); 
        title = "z validty: mean", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    ); push!(ps, p)

    p = heatmap(beta_orders, beta_orders, just_valids(vg_board, vg_valid_board); 
        title = "vg valid: mean", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    ); push!(ps, p)

    p = heatmap(beta_orders, beta_orders, just_valids(H_board, vg_valid_board .* z_valid_board); 
        title = "total valid: H", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    ); push!(ps, p)

    p = heatmap(beta_orders, beta_orders, just_valids(H_board, z_valid_board);
        title = "z validty: H", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    ); push!(ps, p)

    p = heatmap(beta_orders, beta_orders, just_valids(H_board, vg_valid_board); 
        title = "vg validty: H", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    ); push!(ps, p)

    p = heatmap(beta_orders, beta_orders, H_board; 
        title = "total H", xlabel = "vg beta", ylabel = "z beta", 
        common_params...
    ); push!(ps, p)

    layout = (2, 3)
    mysavefig(ps, "feasible_maps"; layout)
    
    # mag vs entropy
    mag(vs) = sqrt(sum(vs .^ 2))

    mags = []
    Hs = []
    for (b1i, beta1) in betas |> enumerate
        for (b2i, beta2) in betas |> enumerate

            z_valid = z_valid_board[b1i, b2i]
            vg_valid = vg_valid_board[b1i, b2i]
            !(z_valid && vg_valid) && continue

            push!(mags, mag([beta1, beta2]))
            push!(Hs, H_board[b1i, b2i])
        end
    end

    p = scatter(mags, Hs; 
        label = "", xlabel = "beta mag", ylabel = "entropy", color = :black
    )
    mysavefig(p, "beta_mag_vs_entropy")

end

## -----------------------------------------------------------------------------------------------
function is_stationary(v, th, w)
    length(v) < w && return false
    @views pv = v[end - w + 1:end]
    m = mean(pv)
    s = std(pv)
    return s <= abs(m) * th
end

## -----------------------------------------------------------------------------------------------
# Find max entropy
ME_GD_DID = (:FIND_MAX_ENTROPY_RESULTS)
let
    M, _ = UJL.load_cache(BEXP_DID)

    biom_idx = M.obj_idx
    LP_cache = Dyn.vgvatp_cache(M)
    z(vatp, vg) = LP_cache[vatp][vg][biom_idx]
    vg(vatp, vg) = vg

    # fake stst params
    cgD_X = 0.3
    M.D = 0.03
    M.X = 1.0
    M.cg = M.X * cgD_X / M.D
    z_avPX = M.D

    # globals
    z_avPME = 0.0
    vg_avPME = 0.0
    PME = Dyn.get_join(M)

    # z_beta, vg_beta = 0.0, 0.0
    # Test
    z_beta, vg_beta = float(rand(-200:200)), float(rand(-200:200))
    z_betas, vg_betas = [z_beta], [vg_beta]

    verb_frec = 10
    gdth = 0.005
    maxiter = 30
    stw = 5
    stth = 0.01

    topiter = 1
    while true

        ## -----------------------------------------------------------------------------------------------
        # find z_beta⁰
        # Gradient descent
        target = z_avPX
        x0 = z_beta
        maxΔx = max(100.0, abs(z_beta) * 0.1)
        x1 = x0 + maxΔx * 0.1

        function z_fun(gdmodel)

            gditer = gdmodel.iter
            z_beta = UJL.gd_value(gdmodel)

            PME = Dyn.get_join!(M, PME) do vatp_, vg_
                exp(z_beta * z(vatp_, vg_) + vg_beta * vg(vatp_, vg_))
            end
            z_avPME = Dyn.ave_over(PME) do vatp, vg
                z(vatp, vg)
            end
            
            err = gdmodel.ϵi
            show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
                gditer == maxiter || err < gdth
            # show_info && begin
            #     @info("Grad Descent ", 
            #         topiter, gditer,  
            #         (z_avPX, z_avPME), 
            #         err, z_beta
            #     ); println()
            # end

            return z_avPME 
        end

        gdmodel = UJL.grad_desc(z_fun; gdth, 
            target, x0, x1, maxΔx, 
            maxiter, 
            verbose = true
        )
        z_beta = UJL.gd_value(gdmodel)
        
        ## -----------------------------------------------------------------------------------------------
        # Check balance at vg_beta = 0.0
        PME = Dyn.get_join!(M, PME) do vatp_, vg_
            vg_beta_ = 0.0
            exp(z_beta * z(vatp_, vg_) + vg_beta_ * vg(vatp_, vg_))
        end
        vg_avPME_at_b0 = Dyn.ave_over(PME) do vatp_, vg_
            vg(vatp_, vg_)
        end
        vg_valid_at_b0 = vg_avPME_at_b0 <= cgD_X

        ## -----------------------------------------------------------------------------------------------
        # if not valid, move vg_beta
        if vg_valid_at_b0
            vg_beta = 0.0
        else
            # Gradient descent
            target = cgD_X * (1.0 - gdth)
            x0 = vg_beta
            maxΔx = max(10.0, abs(vg_beta) * 0.1)
            x1 = x0 + maxΔx * 0.1
            
            function vg_fun(gdmodel)

                gditer = gdmodel.iter
                vg_beta = UJL.gd_value(gdmodel)

                PME = Dyn.get_join!(M, PME) do vatp_, vg_
                    exp(z_beta * z(vatp_, vg_) + vg_beta * vg(vatp_, vg_))
                end
                vg_avPME = Dyn.ave_over(PME) do vatp_, vg_
                    vg(vatp_, vg_)
                end
                
                err = gdmodel.ϵi
                show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
                    gditer == maxiter || err < gdth
                # show_info && begin
                #     @info("Grad Descent ", 
                #         topiter, gditer,  
                #         (cgD_X, vg_avPME), 
                #         err, z_beta
                #     ); println()
                # end

                return vg_avPME 
            end

            gdmodel = UJL.grad_desc(vg_fun; gdth, 
                target, x0, x1, maxΔx, 
                maxiter, 
                verbose = true
            )

            vg_beta = UJL.gd_value(gdmodel)
        end

        push!(z_betas, z_beta)
        push!(vg_betas, vg_beta)

        conv = (vg_valid_at_b0 && gdmodel.ϵi < gdth) || 
            (is_stationary(z_betas, stth, stw) && is_stationary(vg_betas, stth, stw))
        
        @info("Finished Round", 
            topiter, conv,
            (z_beta, vg_beta),
            (vg_avPME, cgD_X),
            (vg_avPME_at_b0, cgD_X),
            (z_avPME, z_avPX),
        ); println()

        conv && break
        topiter += 1

        # Test
        topiter > 20 && break
    end

    # saving results
    dat = (;M, z_betas, vg_betas, cgD_X, z_avPX)
    UJL.save_cache(ME_GD_DID, dat)
end

## -----------------------------------------------------------------------------------------------
# gd plots
let
    M, z_betas, vg_betas = UJL.load_cache(ME_GD_DID)

    z_avPX = M.D
    cgD_X = M.cg * M.D/ M.X

    biom_idx = M.obj_idx
    LP_cache = Dyn.vgvatp_cache(M)
    z(vatp, vg) = LP_cache[vatp][vg][biom_idx]
    vg(vatp, vg) = vg

    common_params = (;
        label = "",
        alpha = 0.7, color = :black,
        xrotation = 35
    )

    ## --------------------------------------------------
    # gd_betas
    ps = Plots.Plot[]
    p = plot(;xlabel = "vg_betas", ylabel = "z_betas")
    scatter!(p, vg_betas, z_betas; 
        common_params...
    )
    plot!(p, vg_betas, z_betas; 
        common_params..., ls = :dash
    ); push!(ps, p)

    ## --------------------------------------------------
    # Entropy vs betas
    mag(vs) = sqrt(sum(vs .^ 2))
    PME = Dyn.get_join(M)
    Hs, mags, z_avPMEs, vg_avPMEs = [], [], [], []
    # colors = palette([:blue, :red], length(z_betas))
    for (z_beta, vg_beta) in zip(z_betas, vg_betas)
        PME = Dyn.get_join!(M, PME) do vatp_, vg_
            exp(z_beta * z(vatp_, vg_) + vg_beta * vg(vatp_, vg_))
        end

        H = entropy(PME)
        m = mag([z_beta, vg_beta])
        push!(Hs, H); push!(mags, m)
        
        z_avPME = Dyn.ave_over(PME) do vatp_, vg_
            z(vatp_, vg_)
        end
        vg_avPME = Dyn.ave_over(PME) do vatp_, vg_
            vg(vatp_, vg_)
        end
        push!(z_avPMEs, z_avPME); push!(vg_avPMEs, vg_avPME)
    end

    # mag vs entropy
    p = plot(;xlabel = "mag", ylabel = "Entropy")
    scatter!(p, mags, Hs; 
        common_params...
    )
    plot!(p, mags, Hs; 
        common_params..., ls = :dash
    ); push!(ps, p)
    
    # mag vs z_avPME
    p = plot(;xlabel = "mag", ylabel = "z_avPME")
    scatter!(p, mags, z_avPMEs; 
        common_params...
    )
    plot!(p, mags, z_avPMEs; 
        common_params..., ls = :dash
    ); push!(ps, p)
    hline!(p, [z_avPX]; common_params...)

    # mag vs vg_avPMEs
    p = plot(;xlabel = "mag", ylabel = "vg_avPME")
    scatter!(p, mags, vg_avPMEs; 
        common_params...
    )
    plot!(p, mags, vg_avPMEs; 
        common_params..., ls = :dash
    ); push!(ps, p)
    hline!(p, [cgD_X]; common_params...)

    ## --------------------------------------------------
    # Walk in beta validity map
    _, beta_orders, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    betas = getbeta.(beta_orders)
    z_th = 0.08
    z_valid_board, vg_valid_board = get_validity_boards(M, betas, z_board, vg_board; z_th)
    # validity maps
    
    function just_valids(dat, validity)
        vals = dat .* validity
        ifelse.(iszero.(vals), NaN, vals)
    end

    valids = just_valids(H_board, z_valid_board .* vg_valid_board)
    p = heatmap(betas, betas, valids; 
        title = "valids: heat H", xlabel = "vg beta", ylabel = "z beta", 
    )

    # gd_betas
    scatter!(p, vg_betas, z_betas; 
        common_params..., 
        color = :red
    )
    scatter!(p, [last(vg_betas)], [last(z_betas)]; 
        label = "", m = (8, :star),
        color = :red
    )

    plot!(p, vg_betas, z_betas; 
        common_params..., ls = :dash
    ); 
    xlim = (-150, 10)
    ylim = (-50, 750)
    plot!(p; xlim, ylim)
    push!(ps, p)
    
    mysavefig(p, "gd_entropy_map_and_path")

    ## --------------------------------------------------
    # saving
    mysavefig(ps, "gd_entropy")
end  