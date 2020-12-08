using Chemostat_Dynamics
using Chemostat_Dynamics.LP_Implement
using ProgressMeter
using Plots
using Serialization
using Base.Threads
using DrWatson

## ---------------------------------------------------------
# SimTools
const CH3_DATA_DIR = joinpath(DATA_DIR, "Ch3Sim")
!isdir(CH3_DATA_DIR) && mkpath(CH3_DATA_DIR)
const CH3_FIGURES_DIR = joinpath(FIGURES_DIR, "Ch3Sim")
!isdir(CH3_FIGURES_DIR) && mkpath(CH3_FIGURES_DIR)

function vrange(L, U, d::Int)
    dx = 10.0^(-d)
    x0 = round(L, RoundUp; digits = d)
    x1 = round(U, RoundDown; digits = d)
    return x0:dx:x1
end

cache_file(name, args...) = joinpath(CACHE_DIR, string(name, hash(args), ".jld"))

function vgvatp_cache(net, θvg, θvatp, 
        vg_idx, vatp_idx, obj_idx; 
        cache_marginf = 0.5)  

    # Open network
    # TODO Detect exchange and open it
    net = deepcopy(net)
    lb, ub = fva(net)
    net.lb .= lb
    net.ub .= ub
    idxs = [vg_idx, vatp_idx]
    margin = cache_marginf .* abs.(net.ub[idxs] .- net.lb[idxs])
    net.lb[idxs] .= net.lb[idxs] .- margin
    net.ub[idxs] .= net.ub[idxs] .+ margin

    # cache = Dict{Tuple{Float64,Float64}, Vector{Float64}}()
    cache = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    # ranges
    vatpL, vatpU = net.lb[vatp_idx], net.ub[vatp_idx]
    vatp_range = vrange(vatpL, vatpU, θvatp)
    vg_ranges = Vector{typeof(vatp_range)}(undef, length(vatp_range))
    for (i, vatp) in enumerate(vatp_range)
        vgL, vgU = fixxing(net, vatp_idx, vatp) do 
            try; return fva(net, vg_idx)
                catch; return 0.0, 0.0
            end
        end
        margin = cache_marginf * abs(vgL - vgU)
        vg_range = vrange(vgL - margin, vgU + margin, θvg)
        vg_ranges[i] = vg_range
    end
    N = sum(length.(values(vg_ranges)))

    # Caching
    prog = Progress(N; desc = "Caching $N flxs ... ", dt = 0.5)
    for (i, vatp) in enumerate(vatp_range)
        lnet = deepcopy(net)
        cache[vatp] = Dict{Float64, Vector{Float64}}()
        for vg in vg_ranges[i]
            fixxing(lnet, vatp_idx, vatp) do 
                fixxing(lnet, vg_idx, vg) do 
                    sol = fba(lnet, obj_idx)
                    cache[vatp][vg] = isempty(sol) ? zeros(length(lnet.c)) : sol            
                end
            end

            next!(prog; showvalues =[
                    (:margin, margin),
                    (:vg, vg),
                    (:vg_range, string(vg_ranges[i], " len: ", length(vg_ranges[i]))),
                    (:vatp, vatp),
                    (:vatp_range, string(vatp_range, " len: ", length(vatp_range))),
                ]
            )   
        end
    end
    finish!(prog)
    return cache
end

## ---------------------------------------------------------
# Sim function
function run_simulation(;
        # Sim params
        θvatp = 2,
        θvg = 3,
        X0 = 1e-3,
        Xmin = 1e-6,
        
        cg = 15.0,
        sg = 15.0,
        Kg = 0.5,
        Vg = 0.5,

        cl = 0.0,
        sl = 0.0,
        Kl = 0.5,
        Vl = 0.0, # 0.1

        D = 1e-2,
        ϵ = 0.0,
        niters = 200,
        damp = 0.9,
        cache_marginf = 0.2
    )

    ## ---------------------------------------------------------
    # model
    net = ToyModel()
    vatp_idx = rxnindex(net, "vatp")
    vg_idx = rxnindex(net, "gt")
    vl_idx = rxnindex(net, "lt")
    obj_idx = rxnindex(net, "biom")
    net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * sg) / (Kg + sg))
    net.ub[vl_idx] = max(net.lb[vl_idx], (Vl * sl) / (Kl + sl))
    
    ## ---------------------------------------------------------
    # Polytope cache
    cfile = cache_file("cache", hash(net), θvatp, θvg, cache_marginf)
    if !isfile(cfile)
        cache = vgvatp_cache(net, θvg, θvatp, 
            vg_idx, vatp_idx, # dimentions
            obj_idx;
            cache_marginf
        )
        serialize(cfile, cache)
    else
        cache = deserialize(cfile)
    end
    
    sl_ts = []
    sg_ts = []
    X_ts = []
    Xb = Dict{Float64, Dict{Float64, Float64}}()

    prog = Progress(niters; desc = "Simulating ... ")
    for it in 1:niters

        ## ---------------------------------------------------------
        # ranges
        vatpL, vatpU = fva(net, vatp_idx)
        vatp_range = vrange(vatpL, vatpU, θvatp)
        vg_ranges = Dict()
        for vatp in vatp_range
            vgL, vgU = fixxing(net, vatp_idx, vatp) do 
                fva(net, vg_idx)
            end
            vg_ranges[vatp] = vrange(vgL, vgU, θvg)
        end

        ## ---------------------------------------------------------
        # X
        if isempty(Xb) 
            for vatp in vatp_range
                Xb[vatp] = Dict{Float64, Float64}() 
                for vg in vg_ranges[vatp]
                    Xb[vatp][vg] = X0
                end
            end
        end

        z = Dict{Float64, Float64}()
        Σvg__X = Dict{Float64, Float64}()
        vgN = Dict{Float64, Int}()

        for vatp in vatp_range
            Σvg__X[vatp] = 0.0
            vg_range = vg_ranges[vatp]
            lXb = get!(Xb, vatp, Dict())
            for vg in vg_range
                lX = get!(lXb, vg, Xmin)
                Σvg__X[vatp] += lX
            end

            vg0 = first(vg_range)
            z[vatp] = cache[vatp][vg0][obj_idx]
            vgN[vatp] = length(vg_range)
        end

        Σ_vatp_vg__z_X = sum(z[vatp]* Σvg__X[vatp] for vatp in vatp_range)
        vatpvgN = sum(values(vgN))
        term2 = ϵ * (Σ_vatp_vg__z_X) / vatpvgN

        for vatp in vatp_range
            term1 = (1 - ϵ) * (z[vatp]* Σvg__X[vatp]) / vgN[vatp]
            lXb = Xb[vatp]
            for vg in vg_ranges[vatp]
                lX = lXb[vg]
                term3 = lX * D
                # update X
                lX = (damp * lX) + (1 - damp) * (lX + term1 + term2 - term3) 
                lXb[vg] = max(Xmin, lX)
            end
        end
        
        ## ---------------------------------------------------------
        # Update Xb
        temp_Xb = Xb
        Xb = Dict{Float64, Dict{Float64, Float64}}()
        for vatp in vatp_range
            lXb = Xb[vatp] = Dict{Float64, Float64}() 
            for vg in vg_ranges[vatp]
                lXb[vg] = temp_Xb[vatp][vg]
            end
        end
        X = sum(sum.(values.(values(Xb))))

        ## ---------------------------------------------------------
        # concs
        Σ_vatp_vg__vg_X = 0.0
        Σ_vatp_vg__vl_X = 0.0
        for vatp in vatp_range
            lcache = cache[vatp]
            lXb = Xb[vatp]
            for vg in vg_ranges[vatp]
                vl = lcache[vg][vl_idx]
                lX = lXb[vg]
                Σ_vatp_vg__vg_X += vg * lX
                Σ_vatp_vg__vl_X += vl * lX
            end
        end
        # update sg
        dsg = -Σ_vatp_vg__vg_X + D*(cg - sg)
        sg = max(0.0, damp * sg + (1 - damp) * (sg + dsg))
        # update sl
        dsl = -Σ_vatp_vg__vl_X + D*(cl - sl)
        sl = max(0.0, damp * sl + (1 - damp) * (sl + dsl))
        

        ## ---------------------------------------------------------
        # Store
        push!(X_ts, X)
        push!(sl_ts, sl)
        push!(sg_ts, sg)

        ## ---------------------------------------------------------
        # Update polytope
        net.ub[vg_idx] = max(net.lb[vg_idx], (Vg * sg) / (Kg + sg))

        ## ---------------------------------------------------------
        # Verbose
        update!(prog, it; showvalues = [
                (:it, it),
                (:X, X),
                (:sg, sg),
                (:dsg, dsg),
                (:sl, sl),
                (:dsl, dsl),
                (:vg_ub, net.ub[vg_idx]),
                (:vatpvgN, vatpvgN),
                (:Xb_len, length(Xb)),
                (:damp, damp)
            ]
        )

    end # for it in 1:niters
    finish!(prog)

    X_ts, sg_ts, sl_ts, Xb
end

## ---------------------------------------------------------

## ---------------------------------------------------------
let
    Ds = collect(10.0.^(-7:0.5:-5))
    for D in Ds
        @info "Doing" D
        P = (;
            θvatp = 2, θvg = 3,
            D, 
            X0 = 1e-5, 
            Xmin = 1e-20, 
            niters = 500, 
            damp = 0.95,
            ϵ = 0.0,
            cache_marginf = 0.1
        )

        X_ts, sg_ts, sl_ts, Xb = run_simulation(;P...)

        fname = savename("Ch3sim_", (;P.D, P.X0, P.niters, P.ϵ))
        simfile = joinpath(CH3_DATA_DIR, string(fname, ".jld"))
        serialize(simfile, (X_ts, sg_ts, sl_ts, Xb, P))
        p1 = plot(xlabel = "time", ylabel = "conc")
        plot!(p1, sg_ts; label = "sg", lw = 3)
        plot!(p1, sl_ts; label = "sl", lw = 3)
        
        p2 = plot(xlabel = "time", ylabel = "X")
        plot!(p2, X_ts; label = "X", lw = 3)
        
        p = plot([p1, p2]...)

        # save fig
        figname = joinpath(CH3_FIGURES_DIR, string(fname, "png"))
        savefig(p, figname)
    end
end

