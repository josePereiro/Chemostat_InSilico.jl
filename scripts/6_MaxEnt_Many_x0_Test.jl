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
fileid = "6"
mysavefig(p, pname; params...) = 
    Dyn.mysavefig(p, pname, Dyn.plotsdir(), fileid; params...)

## -----------------------------------------------------------------------------------------------
# MaxEnt Full Polytope
ME_FULL_POL_DID = (:MAX_ENTROPY_FULL_POLYTOPE_RESULTS)
INDEX_DID = (ME_FULL_POL_DID, :INDEX)
BEXP_DID = (:BETA_SPACE_EXPLORATION_RESULTS)
let
    # error("Comment this line to continue")

    M0, _ = UJL.load_cache(BEXP_DID)

    obj_idx = M0.obj_idx
    LP_cache = Dyn.vgvatp_cache(M0)
    z(vatp, vg) = LP_cache[vatp][vg][obj_idx]
    vg(vatp, vg) = vg

    # fake stst params
    cgD_X = 0.3
    M0.D = 0.03
    M0.X = 1.0
    M0.cg = M0.X * cgD_X / M0.D
    Dyn.check_params(M0)
    z_avPX = M0.D
    
    no_tests = 15
    z_beta0_pool = [float(rand(-200:200, no_tests)); 0.0] |> unique! |> sort!
    vg_beta0_pool = [float(rand(-200:200, no_tests)); 0.0] |> unique! |> sort!
    UJL.save_cache(INDEX_DID, (;z_beta0_pool, vg_beta0_pool))

    to_iter = zip(z_beta0_pool, vg_beta0_pool) |> enumerate |> collect
    @threads for (testi, (z_beta0, vg_beta0)) in to_iter
        
        # globals
        thid = threadid()
        M = deepcopy(M0)
        Dyn.check_params(M0)
        z_avPME = 0.0
        vg_avPME = 0.0

        z_beta, vg_beta = z_beta0, vg_beta0
        z_betas, vg_betas = [z_beta], [vg_beta]
        z_avPMEs, vg_avPMEs = [z_avPME], [vg_avPME]

        PME = Dyn.get_join(M) do vatp_, vg_
            exp(z_beta * z(vatp_, vg_) + vg_beta * vg(vatp_, vg_))
        end
        Hs = [Dyn.entropy(PME)]

        verb_frec = 10
        gdth = 0.005
        maxiter = 100
        stw = 5
        stth = 0.01

        ## -----------------------------------------------------------------------------------------------
        function gd_core_fun(gdmodel; msg)
            PME = Dyn.get_join!(M, PME) do vatp_, vg_
                exp(z_beta * z(vatp_, vg_) + vg_beta * vg(vatp_, vg_))
            end
            z_avPME = Dyn.ave_over(PME) do vatp_, vg_
                z(vatp_, vg_)
            end
            vg_avPME = Dyn.ave_over(PME) do vatp_, vg_
                vg(vatp_, vg_)
            end
            
            err = gdmodel.ϵi
            show_info = gditer == 1 || rem(gditer, verb_frec) == 0 || 
                gditer == maxiter || err < gdth
            show_info && begin
                @info(msg, 
                    topiter, gditer,  
                    (z_avPME, z_avPX), 
                    (vg_avPME, cgD_X), 
                    err, z_beta, vg_beta, thid
                ); println()
            end
        end


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

                msg = "z grad Descent "
                gd_core_fun(gdmodel; msg)

                return z_avPME 
            end

            gdmodel = UJL.grad_desc(z_fun; gdth, 
                target, x0, x1, maxΔx, 
                maxiter, verbose = false
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

                    msg = "vg grad Descent "
                    gd_core_fun(gdmodel; msg)

                    return vg_avPME 
                end

                gdmodel = UJL.grad_desc(vg_fun; gdth, 
                    target, x0, x1, maxΔx, 
                    maxiter, verbose = false
                )

                vg_beta = UJL.gd_value(gdmodel)
            end

            push!(z_betas, z_beta); push!(vg_betas, vg_beta)
            push!(z_avPMEs, z_avPME); push!(vg_avPMEs, vg_avPME)
            push!(Hs, Dyn.entropy(PME))

            conv = (vg_valid_at_b0 && gdmodel.ϵi < gdth) || 
                (UJL.is_stationary(z_betas, stth, stw) && UJL.is_stationary(vg_betas, stth, stw))
            
            @info("Finished Round", 
                topiter, conv,
                (z_beta, vg_beta),
                (vg_avPME, cgD_X),
                (vg_avPME_at_b0, cgD_X),
                (z_avPME, z_avPX),
                thid
            ); println()

            conv && break

            topiter += 1
            topiter > maxiter && break
        end

        # saving results
        did = (ME_FULL_POL_DID, z_beta0, vg_beta0)
        dat = (;M, z_betas, vg_betas, z_avPMEs, vg_avPMEs, Hs, cgD_X, z_avPX)
        UJL.save_cache(did, dat)
        println()

    end # for testi
end

## -----------------------------------------------------------------------------------------------
# Plots

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

## --------------------------------------------------
# gd plots
let
    all_paths_p = nothing
    all_paths_box = Dict(:lx => Inf, :ux => -Inf, :ly => Inf, :uy => -Inf)
    z_beta0_pool, vg_beta0_pool = UJL.load_cache(INDEX_DID)
    for (z_beta0, vg_beta0) in zip(z_beta0_pool, vg_beta0_pool)

        ## --------------------------------------------------
        # load dat
        did = (ME_FULL_POL_DID, z_beta0, vg_beta0)
        dat = UJL.load_cache(did)
        M, z_betas, vg_betas, z_avPMEs, vg_avPMEs, Hs, cgD_X, z_avPX = dat

        ## --------------------------------------------------
        # Plots
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
        mags = map(zip(vg_betas, z_betas)) do (z_beta_, vg_beta_)
            mag([z_beta_, vg_beta_])
        end
        @assert length(mags) == length(Hs)

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
        _, _, betas, H_board, vg_board, z_board = UJL.load_cache(BEXP_DID)

        z_th = 0.08
        z_valid_board, vg_valid_board = get_validity_boards(M, betas, z_board, vg_board; z_th)
        
        # validity maps
        valids = just_valids(H_board, z_valid_board .* vg_valid_board)
        p = heatmap(betas, betas, valids; 
            title = "valids: heat H", xlabel = "vg beta", ylabel = "z beta", 
        )
        isnothing(all_paths_p) && (all_paths_p = deepcopy(p))

        # gd_betas
        for p_ in [p, all_paths_p]
            scatter!(p_, vg_betas, z_betas; 
                common_params..., 
                color = :red
            )
            scatter!(p_, [last(vg_betas)], [last(z_betas)]; 
                label = "", m = (8, :star),
                color = :red
            )
            plot!(p_, vg_betas, z_betas; 
                common_params..., ls = :dash
            ); 

            
        end
        xlim = (min(-150, vg_beta0 - 50), max(10, vg_beta0 + 50))
        ylim = (min(-50, z_beta0 - 50), max(750, z_beta0 + 50))
        plot!(p; xlim, ylim)
        push!(ps, p)
            
        all_paths_box[:lx] = min(all_paths_box[:lx], first(xlim))
        all_paths_box[:ux] = max(all_paths_box[:ux], last(xlim))
        all_paths_box[:ly] = min(all_paths_box[:ly], first(ylim))
        all_paths_box[:uy] = max(all_paths_box[:uy], last(ylim))
        
        mysavefig(p, "gd_entropy_map_and_path"; z_beta0, vg_beta0)
        
        ## --------------------------------------------------
        # saving
        mysavefig(ps, "gd_entropy"; z_beta0, vg_beta0)
        
    end # for (z_beta0, vg_beta0)
    
    xlim = (all_paths_box[:lx] - 50, all_paths_box[:ux] + 50)
    ylim = (all_paths_box[:ly] - 50, all_paths_box[:uy] + 50)
    plot!(all_paths_p; xlim, ylim)
    mysavefig(all_paths_p, "gd_entropy_all_paths")
end  
