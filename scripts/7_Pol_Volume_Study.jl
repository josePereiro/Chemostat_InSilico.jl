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
fileid = "7"
mysavefig(p, pname; params...) = 
    Dyn.mysavefig(p, pname, Dyn.plotsdir(), fileid; params...)

## -----------------------------------------------------------------------------------------------
@time let
    bins = 10
    M = Dyn.SimModel(;fill_Xb = false)
    Ds = range(0.01, 0.048; length = bins)
    cgD_Xs = range(0.01, 0.5; length = bins)
    
    psize = Dyn.poly_vol_board(M, Ds,cgD_Xs; bins)
    
    p = heatmap(Ds, cgD_Xs, psize'; 
        title = "Polytope volume", label = "", 
        xlabel = "D", ylabel = "cgD_X"
    )
    
    mysavefig(p, "plo_volumen"; bins)
    
    # invert for 
    ipsize = similar(psize)
    for (i, _) in enumerate(eachcol(psize))
        ipsize[:, end - (i - 1)] .= psize[:, i]
    end

    yticks = UJL.get_ticks(reverse(cgD_Xs); l = 5) do cgD_X
        -round((maximum(cgD_Xs) - cgD_X); digits = 3)
    end
    
    p = heatmap(Ds, cgD_Xs, ipsize'; 
        title = "Polytope volume", label = "", 
        xlabel = "D", ylabel = "cgD_X", 
        yticks
    )
    
    mysavefig(p, "inv_plo_volumen"; bins)
end

## -----------------------------------------------------------------------------------------------
     