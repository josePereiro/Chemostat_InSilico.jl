import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

import Chemostat_InSilico
import Chemostat_InSilico: CH3_DATA_DIR
const InLP = Chemostat_InSilico.LP_Implement
const InU = Chemostat_InSilico.Utilities

using ProgressMeter
using Plots
using Serialization
using Base.Threads
using Dates
using Base.Threads
using BenchmarkTools
using Statistics

## ---------------------------------------------------------
dat_dir = joinpath(CH3_DATA_DIR, "Ch3___2020-12-20T09:16:19.82")
D = Dict()
for (root, _, files) in walkdir(dat_dir)
    for file in files
        TS, M = deserialize(joinpath(dat_dir, file))
        @show M.D
        D[M.D] = (TS, M)
    end
end

## ---------------------------------------------------------
