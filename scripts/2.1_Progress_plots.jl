import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const Dyn = InCh.Dynamic
    
    # Test
    using Plots
    import GR
    GR.inline("png")

    using Serialization
    using Dates

    import UtilsJL
    const UJL = UtilsJL
    import FileIO
end


## -----------------------------------------------------------------------------------------------
# Plot progress coroutine
const DATA_DIR = InCh.datadir()
const DATA_FILE_PREFFIX = "dyn_dat"
const PROG_FIG_DIR = Dyn.plotsdir("progress")
mkpath(PROG_FIG_DIR)

## -----------------------------------------------------------------------------------------------
# Progress plots
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
            try
                status, TS, M = deserialize(datfile)
                # TODO handle naming with DW
                hist_fname = string(string(fname, "_tslen=",length(TS)), ".png")
                hist_figfile = joinpath(PROG_FIG_DIR, hist_fname)

                @info("Updating progress", fname, now()); println()
                
                p = Dyn.plot_res(M, TS)
                savefig(p, curr_figfile)
                savefig(p, hist_figfile)

            catch err
                @warn("ERROR", err); println()
            end
        end
    end
    sleep(rand(3:8))
end