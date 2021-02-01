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

while true
    
    @info "Scanning " DATA_DIR now()
    println()

    for file in readdir(DATA_DIR)
        !startswith(file, DATA_FILE_PREFFIX) && continue
    
        # dat file
        datfile = joinpath(DATA_DIR, file)

        # fig file
        fname = replace(file, DATA_FILE_PREFFIX => "fig")
        fname, _ = splitext(fname)
        fname = string(fname, ".png")
        figfile = joinpath(PROG_FIG_DIR, fname)

        if !isfile(figfile) || mtime(datfile) > mtime(figfile)
            # save fig
            status, TS, M = deserialize(datfile)
            
            @info "Updating progress" fname now()
            println()
            
            p = InLP.plot_res(M, TS)
            savefig(p, figfile)
        else
            @info "Up to date" fname now()
            println()
        end
    end
    sleep(rand(3:8))
end