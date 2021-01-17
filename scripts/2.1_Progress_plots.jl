import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_InSilico")

@time begin

    import Chemostat_InSilico
    const InCh = Chemostat_InSilico
    const InLP = InCh.LP_Implement
    const InU = InCh.Utilities
    
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
const MTIMES = Dict()

while true
    
    @info "Scanning " DATA_DIR now()
    println()

    for file in readdir(DATA_DIR)
        !startswith(file, DATA_FILE_PREFFIX) && continue
        
        cfile = joinpath(DATA_DIR, file)

        currt = mtime(cfile)
        lastt = get!(MTIMES, file, -Inf)

        # save fig
        if lastt != currt
            status, TS, M = deserialize(cfile)
            fname = replace(file, DATA_FILE_PREFFIX => "fig")
            fname, _ = splitext(fname)
            fname = string(fname, ".png")
            ffile = joinpath(PROG_FIG_DIR, fname)
            
            @info "Plotting progress" fname now()
            println()
            
            p = InLP.plot_res(M, TS)
            savefig(p, ffile)
            MTIMES[file] = currt
        end
    end
    sleep(rand(3:8))
end