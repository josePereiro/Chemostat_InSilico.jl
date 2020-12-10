struct SimModel
    # net
    net::MetNet              # auxiliar metabolic network
    vatp_idx::Int
    vg_idx::Int
    vl_idx::Int
    obj_idx::Int

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

    # Chemostat state
    sl_ts::Vector{Float64}
    sg_ts::Vector{Float64}
    X_ts::Vector{Float64}
    Xb::Dict{Float64, Dict{Float64, Float64}}

    function SimModel(;net = ToyModel(), θvatp = 2, θvg = 3, 
            X0 = 1e-5, Xmin = 1e-20, D = 1e-2, 
            cg = 15.0, sg0 = 15.0, Kg = 0.5, Vg = 0.5, 
            cl = 0.0, sl0 = 0.0, Kl = 0.5, Vl = 0.0, 
            ϵ = 0.0, 
            niters = 200, damp = 0.9,
            sg = Ref{Float64}(sg0), sl = Ref{Float64}(sl0),
            sl_ts = [], sg_ts = [], X_ts = [],
            Xb = Dict{Float64, Dict{Float64, Float64}}(),
            vatp_ider = "vatp",
            vg_ider = "gt",
            vl_ider = "lt",
            obj_ider = "biom"
        ) 

        vatp_idx = rxnindex(net, vatp_ider)
        vg_idx = rxnindex(net, vg_ider)
        vl_idx = rxnindex(net, vl_ider)
        obj_idx = rxnindex(net, obj_ider)

        new(net, vatp_idx, vg_idx, vl_idx, obj_idx,
            θvatp, θvg, X0, Xmin, D, 
            cg, sg0, Kg, Vg, 
            cl, sl0, Kl, Vl, 
            ϵ, niters, damp,  
            sl_ts, sg_ts, X_ts, Xb
        )
    end
end

# function SimModel(M::SimModel; kwargs...)
#     d = Dict{Symbol, Any}()
#     for k in fieldnames(SimModel)
#         d[k] = getfield(M, k)
#     end
#     for (k, v) in kwargs
#         d[k] = v
#     end
#     # return SimModel(;d...)
#     d
# end 

## ---------------------------------------------------------  
function vrange(L, U, d::Int)
    dx = 10.0^(-d)
    x0 = round(L, RoundUp; digits = d)
    x1 = round(U, RoundDown; digits = d)
    return x0:dx:x1
end

function vatpvg_ranges(M::SimModel; 
        vatp_margin::Float64 = 0.0,
        vg_margin::Float64 = 0.0,
    )
    # ranges
    vatpL, vatpU = fva(M.net, M.vatp_idx)
    vatp_range = vrange(vatpL - vatp_margin, vatpU + vg_margin, M.θvatp)
    vg_ranges = Vector{typeof(vatp_range)}(undef, length(vatp_range))
    for (i, vatp) in enumerate(vatp_range)
        vgL, vgU = fixxing(M.net, M.vatp_idx, vatp) do 
            fva(M.net, M.vg_idx)
        end
        @inbounds vg_ranges[i] = vrange(vgL - vg_margin, vgU + vg_margin, M.θvg)
    end
    return vatp_range, vg_ranges
end

## ---------------------------------------------------------  
cache_file(M::SimModel) = joinpath(CACHE_DIR, 
    mysavename("cache_", "jld"; M.θvatp, M.θvg))

function get_cache(M::SimModel)
    cfile = cache_file(M)
    if isfile(cfile)
        cache = deserialize(cfile)
    else
        cache = vgvatp_cache(M)
        serialize(cfile, cache)
    end
    return cache
end
    

## ---------------------------------------------------------  
function vgvatp_cache(M::SimModel)

    # Open network
    M = deepcopy(M)
    net = M.net
    
    # ranges
    vatp_range, vg_ranges = vatpvg_ranges(M::SimModel; 
        vatp_margin = 10.0^(-M.θvatp),
        vg_margin = 10.0^(-M.θvg)
    )
    N = sum(length.(values(vg_ranges)))

    # Caching
    cache = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    prog = Progress(N; desc = "Caching (N = $N)  ... ", dt = 0.5)
    c = 0
    for (i, vatp) in enumerate(vatp_range)
        lnet = deepcopy(net)
        cache[vatp] = Dict{Float64, Vector{Float64}}()
        @inbounds vg_range = vg_ranges[i]
        for vg in vg_range
            fixxing(lnet, M.vatp_idx, vatp) do 
                fixxing(lnet, M.vg_idx, vg) do 
                    sol = fba(lnet, M.obj_idx)
                    cache[vatp][vg] = sol            
                end
            end
            c += 1
        end
        update!(prog, c; showvalues = [
                    (:vatp_range, string(vatp_range, " len: ", length(vatp_range))),
                    (:vg_range, string(vg_range, " len: ", length(vg_range))),
                ]
            )   
    end
    finish!(prog)
    return cache
end

# ---------------------------------------------------------
# Sim function
function run_simulation!(M::SimModel; 
        at_iter = (it, P, O) -> P,
        cache = get_cache(M),
        verbose = true
    )

    ## ---------------------------------------------------------
    # prepare model
    net = M.net
    vatp_idx, vg_idx, vl_idx, obj_idx = M.vatp_idx, M.vg_idx, M.vl_idx, M.obj_idx
    net.ub[vg_idx] = max(net.lb[vg_idx], (M.Vg * M.sg0) / (M.Kg + M.sg0))
    net.ub[vl_idx] = max(net.lb[vl_idx], (M.Vl * M.sl0) / (M.Kl + M.sl0))

    # cache
    if isnothing(cache)
        cache = vgvatp_cache(M)
    end

    # Out
    sl_ts, sg_ts, X_ts, Xb = M.sl_ts, M.sg_ts, M.X_ts, M.Xb
    sg, sl = M.sg0, M.sl0
    
    verbose && (prog = Progress(M.niters; desc = "Simulating ... "))
    for it in 1:M.niters

        ## ---------------------------------------------------------
        # ranges
        vatp_range, vg_ranges = vatpvg_ranges(M::SimModel)

        ## ---------------------------------------------------------
        # X
        if isempty(Xb) 
            for (i, vatp) in enumerate(vatp_range)
                Xb[vatp] = Dict{Float64, Float64}() 
                for vg in @inbounds vg_ranges[i]
                    Xb[vatp][vg] = M.X0
                end
            end
        end

        z = Dict{Float64, Float64}()
        Σvg__X = Dict{Float64, Float64}()
        vgN = Dict{Float64, Int}()

        for (i, vatp) in enumerate(vatp_range)
            Σvg__X[vatp] = 0.0
            @inbounds vg_range = vg_ranges[i]
            lXb = get!(Xb, vatp, Dict())
            for vg in vg_range
                lX = get!(lXb, vg, M.Xmin)
                Σvg__X[vatp] += lX
            end

            vg0 = first(vg_range)
            z[vatp] = cache[vatp][vg0][obj_idx]
            vgN[vatp] = length(vg_range)
        end

        Σ_vatp_vg__z_X = sum(z[vatp]* Σvg__X[vatp] for vatp in vatp_range)
        vatpvgN = sum(values(vgN))
        term2 = M.ϵ * (Σ_vatp_vg__z_X) / vatpvgN

        for (i, vatp) in enumerate(vatp_range)
            term1 = (1 - M.ϵ) * (z[vatp]* Σvg__X[vatp]) / vgN[vatp]
            lXb = Xb[vatp]
            for vg in @inbounds vg_ranges[i]
                lX = lXb[vg]
                term3 = lX * M.D
                # update X
                lX = (M.damp * lX) + (1 - M.damp) * (lX + term1 + term2 - term3) 
                lXb[vg] = max(M.Xmin, lX)
            end
        end
        
        ## ---------------------------------------------------------
        # Update Xb (let unfeasible out)
        temp_Xb = deepcopy(Xb)
        empty!(Xb)
        for (i, vatp) in enumerate(vatp_range)
            lXb = Xb[vatp] = Dict{Float64, Float64}() 
            for vg in @inbounds vg_ranges[i]
                lXb[vg] = temp_Xb[vatp][vg]
            end
        end
        X = sum(sum.(values.(values(Xb))))

        ## ---------------------------------------------------------
        # concs
        Σ_vatp_vg__vg_X = 0.0
        Σ_vatp_vg__vl_X = 0.0
        for (i, vatp) in enumerate(vatp_range)
            lcache = cache[vatp]
            lXb = Xb[vatp]
            for vg in @inbounds vg_ranges[i]
                vl = lcache[vg][vl_idx]
                lX = lXb[vg]
                Σ_vatp_vg__vg_X += vg * lX
                Σ_vatp_vg__vl_X += vl * lX
            end
        end
        # update sg
        dsg = -Σ_vatp_vg__vg_X + M.D * (M.cg - sg)
        sg = max(0.0, M.damp * sg + (1 - M.damp) * (sg + dsg))
        # update sl
        dsl = -Σ_vatp_vg__vl_X + M.D * (M.cl - sl)
        sl = max(0.0, M.damp * sl + (1 - M.damp) * (sl + dsl))
        

        ## ---------------------------------------------------------
        # Store
        push!(X_ts, X)
        push!(sl_ts, sl)
        push!(sg_ts, sg)

        ## ---------------------------------------------------------
        # Update polytope
        net.ub[vg_idx] = max(net.lb[vg_idx], (M.Vg * sg) / (M.Kg + sg))

        ## ---------------------------------------------------------
        # Verbose
        verbose && update!(prog, it; showvalues = [
                (:it, it),
                (:D, M.D),
                (:damp, M.damp),
                (": ------ ", "------------------"),
                (:X, X),
                (:sg, sg),
                (:dsg, dsg),
                (:sl, sl),
                (:dsl, dsl),
                (:vg_ub, net.ub[vg_idx]),
                (:vatpvgN, vatpvgN),
                (:Xb_len, sum(length.(values(Xb)))),
            ]
        )

    end # for it in 1:niters
    verbose && finish!(prog)

    return M
end