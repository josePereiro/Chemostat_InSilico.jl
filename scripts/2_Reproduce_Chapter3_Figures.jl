using Chemostat_Dynamics
using Chemostat_Dynamics.LP_Implement
using ProgressMeter
using Plots
using Serialization
using Base.Threads
using DrWatson
using Printf

## ---------------------------------------------------------
# Sim Tools
sci(n) = @sprintf("%0.1e", n)
function mysavename(name, ext = ""; c...)
    d = Dict{Symbol, Any}(c)
    for (k, v) in d
        if v isa AbstractFloat
            d[k] = sci(v)
        end
    end
    savename(name, d, ext)
end

## ---------------------------------------------------------
# SimTools
let
    lagtime = 1000
    etime = 5000
    Ds = collect(10.0.^(-7:0.5:-2))
    bkps = collect(lagtime .+ etime * (1:length(Ds)))
    bkpi = 1

    P = SimParams(;
        θvatp = 2, θvg = 3,
        D = Ds[bkpi], 
        X0 = 1e-6, 
        Xmin = 1e-20, 
        niters = lagtime + etime * length(Ds) + etime,
        damp = 0.95,
        ϵ = 0.0,
        sg0 = 5.0,
        cache_marginf = 0.1
    )

    function at_iter(it, P) 
        next_bkp = bkps[bkpi]
        (it != next_bkp || bkpi == length(Ds)) && return P
        bkpi += 1
        new_D = Ds[bkpi]
        return SimParams(P; D = new_D)
    end
    X_ts, sg_ts, sl_ts, Xb = run_simulation(P; at_iter)

    fname = mysavename("Ch3sim_"; P.X0, P.niters, P.ϵ)
    @show fname
    simfile = joinpath(CH3_DATA_DIR, string(fname, ".jld"))
    serialize(simfile, (X_ts, sg_ts, sl_ts, Xb, P, Ds))
    # X_ts, sg_ts, sl_ts, Xb, P, Ds = deserialize(simfile)
    p1 = plot(xlabel = "time", ylabel = "conc")
    plot!(p1, log10.(sg_ts); label = "sg", lw = 3)
    plot!(p1, log10.(sl_ts); label = "sl", lw = 3)
    
    p2 = plot(xlabel = "time", ylabel = "X")
    plot!(p2, log10.(X_ts); label = "X", lw = 3)
    
    p = plot([p1, p2]...)

    # save fig
    figname = joinpath(CH3_FIGURES_DIR, string(fname, ".png"))
    savefig(p, figname)
end
