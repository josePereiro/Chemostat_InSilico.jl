import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement
    
    # Test
    using Plots
    import GR
    GR.inline("png")

    using Serialization
    using Dates

end

## -----------------------------------------------------------------------------------------------
# Plot progress coroutine
const DATA_DIR = InCh.DYN_DATA_DIR
const DATA_FILE_PREFFIX = "dyn_dat"
const PROG_FIG_DIR = joinpath(InCh.DYN_FIGURES_DIR, "progress")

## -----------------------------------------------------------------------------------------------
using Random
using ProgressMeter


# ## -----------------------------------------------------------------------------------------------
# let
    
#     # dat file
#     datname = "dyn_dat_D=0.01_Vl=0.0e+00_cg=15_Δt=0.5_σ=0.01_τ=0.0e+00_ϵ=0.5.jls"
#     datfile = joinpath(DATA_DIR, datname)

#     # fig file
#     !isdir(PROG_FIG_DIR) && mkpath(PROG_FIG_DIR)
#     figfile = joinpath(PROG_FIG_DIR, "test.png")

#     status, TS, M = deserialize(datfile)
#     global M0 = M

#     p =  plot(xlabel = "vatp", ylabel = "vg")
#     InLP.plot_polborder!(p, M)
#     plot_poldist!(p, M; 
#         static_th = 0.2, 
#         pkwargs = Dict(), 
#         skwargs = Dict(:alpha => 1.0)
#     )
#     savefig(p, figfile)
    
# end
# -----------------------------------------------------------------------------------------------
while true
    
    @info("Scanning ", DATA_DIR, now()); println()

    for file in readdir(DATA_DIR)
        !startswith(file, DATA_FILE_PREFFIX) && continue
    
        # dat file
        datfile = joinpath(DATA_DIR, file)

        # fig file
        !isdir(PROG_FIG_DIR) && mkpath(PROG_FIG_DIR)
        fname = replace(file, DATA_FILE_PREFFIX => "fig")
        fname, _ = splitext(fname)
        curr_fname = string(fname, ".png")
        curr_figfile = joinpath(PROG_FIG_DIR, curr_fname)
        
        if !isfile(curr_figfile) || mtime(datfile) > mtime(curr_figfile)
            # save fig
            try
                status, TS, M = deserialize(datfile)
                hist_fname = string(string(fname, "_tslen_",length(TS)), ".png")
                hist_figfile = joinpath(PROG_FIG_DIR, hist_fname)

                @info("Updating progress", fname, now()); println()
                
                p = InLP.plot_res(M, TS)
                savefig(p, hist_fname)
                savefig(p, hist_figfile)
            catch err
                @warn("ERROR", err); println()
            end
        end
    end
    sleep(rand(3:8))
end