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

getbeta(order) = sign(order) * 10.0 ^ abs(order)

## -----------------------------------------------------------------------------------------------
# Compute beta maps
BEXP_DID = (:BETA_SPACE_EXPLORATION_RESULTS)
DTASK = @async let
    error("Comment this line to continue")

    # base model
    M = Dyn.SimModel()

    obj_idx = M.obj_idx
    LP_cache = Dyn.vgvatp_cache(M)
    z(vatp, vg) = LP_cache[vatp][vg][obj_idx]
    vg(vatp, vg) = vg

    beta_orders = -3:0.1:3
    betas = getbeta.(beta_orders)

    H_board = zeros(length(beta_orders), length(beta_orders))
    vg_board = similar(H_board)
    z_board = similar(H_board)
    prog = Progress(length(H_board))

    Ch = Channel(nthreads()) do Ch_
        for (z_bi, z_beta) in enumerate(betas)
            for (vg_bi, vg_beta) in enumerate(betas)
                put!(Ch_, (z_bi, z_beta, vg_bi, vg_beta))
            end
        end
    end
    
    @threads for _ in 1:nthreads()
        thid = threadid()
        isinteractive() && thid == 1 && continue # Let main free
        
        local M0 = deepcopy(M)
        local PME = Dyn.get_join(M0)
        for (z_bi, z_beta, vg_bi, vg_beta) in Ch
            
            try;
                PME = Dyn.get_join!(M0, PME) do _vatp, _vg
                    exp(z_beta * z(_vatp, _vg) + vg_beta * vg(_vatp, _vg))
                end
                H = UJL.entropy(PME)

                H_board[z_bi, vg_bi] = H
                vg_board[z_bi, vg_bi] = Dyn.ave_over(PME) do vatp, vg
                    vg
                end
                z_board[z_bi, vg_bi] = Dyn.ave_over(PME) do vatp, vg
                    z(vatp, vg)
                end

                next!(prog; showvalues = [
                        (:z_beta, z_beta),    
                        (:vg_beta, vg_beta),    
                        (:H, H),    
                        (:thid, thid),    
                    ]
                )
            catch err
                @warn("Error at", z_beta, vg_beta)
                rethrow(err)
            end
        end

    end # for thid

    # save
    dat = (;M, beta_orders, betas,  H_board, vg_board, z_board)
    UJL.save_cache(BEXP_DID, dat; verbose = true)

end 

## -----------------------------------------------------------------------------------------------
# raw plots
let
    M, beta_orders, betas, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)

    ticks = UJL.get_ticks(beta_orders; l = 11) do order
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
    mysavefig(ps, "raw-maps")
end

## -----------------------------------------------------------------------------------------------
# entropy plots
let
    M, beta_orders, betas, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    maxbeta = maximum(betas)

    xticks = UJL.get_ticks(eachindex(beta_orders); l = 11) do idx
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

    mysavefig(ps, "entropy-vs-beta")
end


## -----------------------------------------------------------------------------------------------
# beta mag vs entropy
let
    M, beta_orders, betas, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)

    mag(vs) = sqrt(sum(vs .^ 2))

    mags = []
    Hs = []
    for (z_bi, z_beta) in betas |> enumerate
        for (vg_bi, vg_beta) in betas |> enumerate
            push!(mags, mag([z_beta, vg_beta]))
            push!(Hs, H_board[z_bi, vg_bi])
        end
    end

    p = scatter(mags, Hs; 
        label = "", xlabel = "beta mag", ylabel = "entropy", color = :black
    )
    mysavefig(p, "beta-mag-vs-entropy")
end

## -----------------------------------------------------------------------------------------------
function get_validity_boards(M, betas, z_board, vg_board; z_th = 0.05)
    cgD_X = M.cg * M.D / M.X
    z_valid_board = similar(z_board, Bool)
    vg_valid_board = similar(z_board, Bool)
    z_th = M.D * z_th
    for (z_bi, z_beta) in betas |> enumerate
        for (vg_bi, vg_beta) in betas |> enumerate
            z = z_board[z_bi, vg_bi]
            vg = vg_board[z_bi, vg_bi]

            # check validity
            z_valid_board[z_bi, vg_bi] = abs(z - M.D) < z_th
            vg_valid_board[z_bi, vg_bi] = vg < cgD_X
        end
    end
    z_valid_board, vg_valid_board
end

function just_valids(dat, validity)
    vals = dat .* validity
    ifelse.(iszero.(vals), NaN, vals)
end

## -----------------------------------------------------------------------------------------------
let
    M, beta_orders, betas, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    betas = getbeta.(beta_orders)

    ticks = UJL.get_ticks(beta_orders; l = 11) do order
        round(getbeta(order))
    end

    common_params = (;
        xticks = ticks,
        yticks = ticks,
        xrotation = 35,
        grid = false,
    )

    # fake stst params
    cgD_Xs = 0.05:0.01:0.35
    Ds = 0.005:0.001:0.04
    Xs = [1.0]
    
    ps = Plots.Plot[]
    for (cgD_X, D, X) in Iterators.product(cgD_Xs, Ds, Xs)
        M.D, M.X, M.cg = D, X, X * cgD_X / D
        
        z_th = 0.08
        z_valid_board, vg_valid_board = get_validity_boards(M, betas, z_board, vg_board; z_th)
        
        @info("Doing: ", (D, X, cgD_X), z_th)
        
        title = string("valid region: H\n", 
            "cg: ", UJL.sci(M.cg), 
            ", D: ", UJL.sci(M.D), 
            ", X: ", UJL.sci.(M.X), 
        )

        p = heatmap(beta_orders, beta_orders, just_valids(H_board, vg_valid_board .* z_valid_board); 
            title, xlabel = "vg beta", ylabel = "z beta", 
            common_params...
        ); 
        # mysavefig(p, "total-valid-H"; M.D, M.X, M.cg)
        push!(ps, p)
    end

    file = Dyn.plotsdir("total-valid-H.gif")
    @time UJL.save_gif(ps, file; fps = 10.0)
end


## -----------------------------------------------------------------------------------------------
# validity plots
let
    M, beta_orders, betas, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)
    betas = getbeta.(beta_orders)

    # fake stst params
    cgD_X = 0.3
    M.D = 0.03
    M.X = 1.0
    M.cg = M.X * cgD_X / M.D
    
    z_th = 0.08
    z_valid_board, vg_valid_board = get_validity_boards(M, betas, z_board, vg_board; z_th)

    ticks = UJL.get_ticks(beta_orders; l = 11) do order
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
    ); 
    mysavefig(p, "total_valid_H"; M.sg, M.D, M.X)
    push!(ps, p)

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
    ); 
    push!(ps, p)

    layout = (2, 3)
    mysavefig(ps, "feasible_maps"; layout)
    
    # mag vs entropy
    mag(vs) = sqrt(sum(vs .^ 2))

    mags = []
    Hs = []
    for (z_bi, z_beta) in betas |> enumerate
        for (vg_bi, vg_beta) in betas |> enumerate

            z_valid = z_valid_board[z_bi, vg_bi]
            vg_valid = vg_valid_board[z_bi, vg_bi]
            !(z_valid && vg_valid) && continue

            push!(mags, mag([z_beta, vg_beta]))
            push!(Hs, H_board[z_bi, vg_bi])
        end
    end

    p = scatter(mags, Hs; 
        label = "", xlabel = "beta mag", ylabel = "entropy", color = :black
    )
    mysavefig(p, "beta_mag_vs_entropy")

end

