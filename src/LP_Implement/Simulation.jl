struct SimParams
    # Sim params
    θvatp::Int               # exactness of vatp discretization dvatp = 10.0^-(θvatp)
    θvg::Int                 # exactness of vatp discretization dvatp = 10.0^-(θvatp)
    X0::Float64              # initial cell density value of each quanta of the polytope
    Xmin::Float64            # minimal allowed cell density value of each quanta of the polytope
    D::Float64               # Chemostat dilution rate
    
    cg::Float64              # feed G concentration
    sg0::Float64             # initial medium G concenration
    Kg::Float64              # g MM constant
    Vg::Float64              # upper g max bound

    cl::Float64              # feed G concentration
    sl0::Float64             # initial medium G concenration
    Kl::Float64              # g MM constant
    Vl::Float64 # 0.1        # upper g max bound

    ϵ::Float64               # mutation rate (%)

    niters::Int              # simulation iteration number
    damp::Float64            # numeric damp
    cache_marginf::Float64   # cache offset

    SimParams(;θvatp = 2, θvg = 3, 
            X0 = 1e-3, Xmin = 1e-20, D = 1e-2, 
            cg = 15.0, sg0 = 15.0, Kg = 0.5, Vg = 0.5, 
            cl = 0.0, sl0 = 0.0, Kl = 0.5, Vl = 0.0, 
            ϵ = 0.0, 
            niters = 200, damp = 0.9, cache_marginf = 0.2
        ) = new(θvatp, θvg, X0, Xmin, D, cg, sg0, Kg, Vg, cl, sl0, Kl, Vl, 
                ϵ, niters, damp, cache_marginf)
end

function SimParams(P::SimParams; kwargs...)
    d = Dict{Symbol, Any}()
    for k in fieldnames(SimParams)
        d[k] = getfield(P, k)
    end
    for (k, v) in kwargs
        d[k] = v
    end
    return SimParams(;d...)
end 

## ---------------------------------------------------------
# TODO: Finish this  
# struct SimOut
#     sl_ts = []
#     sg_ts = []
#     X_ts = []
#     Xb = Dict{Float64, Dict{Float64, Float64}}()
# end
## ---------------------------------------------------------  
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

        end
        next!(prog; showvalues = [
                    (:margin, margin),
                    (:vg_range, string(vg_ranges[i], " len: ", length(vg_ranges[i]))),
                    (:vatp_range, string(vatp_range, " len: ", length(vatp_range))),
                ]
            )   
    end
    finish!(prog)
    return cache
end

## ---------------------------------------------------------
# Sim function
function run_simulation(P::SimParams; at_iter = (it, P) -> P)

    ## ---------------------------------------------------------
    # model
    net = ToyModel()
    vatp_idx = rxnindex(net, "vatp")
    vg_idx = rxnindex(net, "gt")
    vl_idx = rxnindex(net, "lt")
    obj_idx = rxnindex(net, "biom")
    net.ub[vg_idx] = max(net.lb[vg_idx], (P.Vg * P.sg0) / (P.Kg + P.sg0))
    net.ub[vl_idx] = max(net.lb[vl_idx], (P.Vl * P.sl0) / (P.Kl + P.sl0))

    ## ---------------------------------------------------------
    # Polytope cache
    cfile = cache_file("cache", hash(net), P.θvatp, P.θvg, P.cache_marginf)
    if !isfile(cfile)
        cache = vgvatp_cache(net, P.θvg, P.θvatp, 
            vg_idx, vatp_idx, # dimentions
            obj_idx;
            P.cache_marginf
        )
        serialize(cfile, cache)
    else
        cache = deserialize(cfile)
    end

    sg = P.sg0
    sl = P.sl0
    # Out
    sl_ts = []
    sg_ts = []
    X_ts = []
    Xb = Dict{Float64, Dict{Float64, Float64}}()

    prog = Progress(P.niters; desc = "Simulating ... ")
    for it in 1:P.niters

        P = at_iter(it, P)::SimParams

        ## ---------------------------------------------------------
        # ranges
        vatpL, vatpU = fva(net, vatp_idx)
        vatp_range = vrange(vatpL, vatpU, P.θvatp)
        vg_ranges = Dict()
        for vatp in vatp_range
            vgL, vgU = fixxing(net, vatp_idx, vatp) do 
                fva(net, vg_idx)
            end
            vg_ranges[vatp] = vrange(vgL, vgU, P.θvg)
        end

        ## ---------------------------------------------------------
        # X
        if isempty(Xb) 
            for vatp in vatp_range
                Xb[vatp] = Dict{Float64, Float64}() 
                for vg in vg_ranges[vatp]
                    Xb[vatp][vg] = P.X0
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
                lX = get!(lXb, vg, P.Xmin)
                Σvg__X[vatp] += lX
            end

            vg0 = first(vg_range)
            z[vatp] = cache[vatp][vg0][obj_idx]
            vgN[vatp] = length(vg_range)
        end

        Σ_vatp_vg__z_X = sum(z[vatp]* Σvg__X[vatp] for vatp in vatp_range)
        vatpvgN = sum(values(vgN))
        term2 = P.ϵ * (Σ_vatp_vg__z_X) / vatpvgN

        for vatp in vatp_range
            term1 = (1 - P.ϵ) * (z[vatp]* Σvg__X[vatp]) / vgN[vatp]
            lXb = Xb[vatp]
            for vg in vg_ranges[vatp]
                lX = lXb[vg]
                term3 = lX * P.D
                # update X
                lX = (P.damp * lX) + (1 - P.damp) * (lX + term1 + term2 - term3) 
                lXb[vg] = max(P.Xmin, lX)
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
        dsg = -Σ_vatp_vg__vg_X + P.D * (P.cg - sg)
        sg = max(0.0, P.damp * sg + (1 - P.damp) * (sg + dsg))
        # update sl
        dsl = -Σ_vatp_vg__vl_X + P.D * (P.cl - sl)
        sl = max(0.0, P.damp * sl + (1 - P.damp) * (sl + dsl))
        

        ## ---------------------------------------------------------
        # Store
        push!(X_ts, X)
        push!(sl_ts, sl)
        push!(sg_ts, sg)

        ## ---------------------------------------------------------
        # Update polytope
        net.ub[vg_idx] = max(net.lb[vg_idx], (P.Vg * sg) / (P.Kg + sg))

        ## ---------------------------------------------------------
        # Verbose
        update!(prog, it; showvalues = [
                (:it, it),
                (:X, X),
                (:D, P.D),
                (:sg, sg),
                (:dsg, dsg),
                (:sl, sl),
                (:dsl, dsl),
                (:vg_ub, net.ub[vg_idx]),
                (:vatpvgN, vatpvgN),
                (:Xb_len, sum(length.(values(Xb)))),
                (:damp, P.damp)
            ]
        )

    end # for it in 1:niters
    finish!(prog)

    X_ts, sg_ts, sl_ts, Xb
end