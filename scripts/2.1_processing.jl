import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

import Chemostat_InSilico
const ChIS = Chemostat_InSilico
const InLP = Chemostat_InSilico.LP_Implement

using ProgressMeter
using Plots

# import GR
# GR.inline("png")

using UtilsJL
using Serialization
using Statistics

## ----------------------------------------------------------------------------
dat_dir = joinpath(ChIS.CH3_DATA_DIR, "LP_Implementation__2020-12-21T18:15:16.132")
@assert isdir(dat_dir)

DAT = DictTree()
for (root, _, files) in walkdir(dat_dir)
    for file in files
        TS, M = deserialize(joinpath(dat_dir, file))
        @info "Loaded" M.D M.ϵ file
        DAT[:TS, M.D, M.ϵ] = TS
        DAT[:M, M.D, M.ϵ] = M
        push!(get!(DAT, [], :Ds), M.D)
        push!(get!(DAT, [], :ϵs), M.ϵ)
    end
    sort!(unique!(DAT[:Ds]))
    sort!(unique!(DAT[:ϵs]))
end

## ----------------------------------------------------------------------------
let
    M = DAT[:M, first(DAT[:Ds]), first(DAT[:ϵs])]
    InLP.vgvatp_cache(M)
end
## ----------------------------------------------------------------------------
# X vs time vs D
let
    f = identity
    marginf = 0.2
    ps = map(DAT[:Ds][1:6]) do D0
    
        Xss = getfield.(DAT[:TS, D0, DAT[:ϵs]], :X_ts) 
        ylim = InLP.lims(marginf, Xss...)
        
        p = plot(xlabel = "time", ylabel = "X", title = string("D: ", UtilsJL.sci(D0)))
        for (ϵ, Xs) in zip(DAT[:ϵs], Xss)
            plot!(p, f.(Xs); 
                label = ϵ, lw = 4, alpha = 0.7, color = Gray(ϵ * 0.8))
        end
        p
    end
    plot(ps...; layout = length(ps), legend = false, titlefont = 12)
end

## ----------------------------------------------------------------------------
using ColorSchemes
let
    plot(rand(10); color = gray(0.4), lw = 6)
end
## ----------------------------------------------------------------------------
# 
let
    p = plot(xlabel = "D", ylabel = "X", title = "Steady State")
    for ϵ in DAT[:ϵs]
        Xs = [DAT[D, ϵ].M.X for D in DAT[:Ds]]
        plot!(p, DAT[:Ds], log10.(Xs .+ 1e-8); label = ϵ, lw = 3, alpha = 0.8)
    end
    p
end

## ----------------------------------------------------------------------------

