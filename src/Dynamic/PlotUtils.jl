## ---------------------------------------------------------
function plot_frees(net::MetNet, frees; 
        δ = 40,
        verbose = true, 
        m = (8, :square),
        pltkwargs...
    )
    freeis = rxnindex(net, frees)
    frees = rxns(net, freeis)

    box = BoxGrid(net, freeis, δ)
    V = length(box)

    # Start collecting, TODO: make threaded
    tolf = 0.01
    tol = Δv(net, freeis) .* tolf

    boxc = 0
    polc = 0

    cid = ("plot_frees", hash(net), δ, hash(freeis))
    polV = Ass.lcache(cid; verbose = false) do
        
        polV_ = [Float64[] for id in frees]

        prog = Progress(V; dt = 0.1, desc = "Collecting... ")
        for freev in box
            fullv = fixxing(net, freeis, freev; tol) do
                fba(net)
            end
            
            # info
            boxc += 1
            verbose && next!(prog; showvalues = [(:polfrac, polc/boxc)])

            # feasible
            isfea = any(fullv .!= 0.0)
            !isfea && continue

            # collect
            for (v, id) in zip(freev, eachindex(frees))
                push!(polV_[id], v)
            end
            
            # info
            polc += 1
            
        end
        verbose && finish!(prog)

        return polV_
    end

    # --------------------------------------------
    println()
    verbose && @info("Ploting")
    ps = Plots.Plot[]
    defaults = (;label = "", m, alpha = 0.5, stroke = false)

    for i in 1:length(frees), j in (i + 1):length(frees)
        freex, freey = frees[i], frees[j]
        polVx, polVy = polV[i], polV[j]
        
        verbose && @show freex, freey
        
        p = plot(;xlabel = freex, ylabel = freey)

        # reduce redundancy
        pool = Set{Tuple{Float64, Float64}}()
        curr_vx = nothing
        for si in sortperm(polVx)
            vx, vy = polVx[si], polVy[si]

            if vx != curr_vx # new vx
                # plot if not empty
                !isempty(pool) && scatter!(p, first.(pool), last.(pool); color = :red,
                    defaults..., pltkwargs...)

                empty!(pool)
                curr_vx = vx
            end
            push!(pool, (vx, vy))
        end
        
        push!(ps, p)
    end
    return ps

end