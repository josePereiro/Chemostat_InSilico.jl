## ---------------------------------------------------------------
let
    iter = enumerate(Iterators.product(Ds, ϵs, cgs))
    Ds_idxs = Dict(D => i for (i, D) in enumerate(Ds))
    ϵs_idxs = Dict(ϵ => i for (i, ϵ) in enumerate(ϵs))

    for cg in cgs
        mat = zeros(length(Ds_idxs), length(ϵs_idxs))
        
        for (D, ϵ) in Iterators.product(Ds, ϵs)
            simparams = (;D, ϵ, cg)
            
            status = get_status(simid, (;D, ϵ, cg))
            (status == UNDONE_SIM_STATUS) && continue

            if status == STST_SIM_STATUS
                S = load_simdat(simid, simparams)
                X = S.X
            elseif status == EXPLODED_SIM_STATUS
                X = 1e4
            elseif status == DEAD_SIM_STATUS
                X = 1e-4
            else
                error("Unexpected status: ", status)
            end

            mat[Ds_idxs[D], ϵs_idxs[ϵ]] = X
        end

        @info("At", cg)
        p = plot(;xlabel = "D", ylabel = "ϵ", title = string("cg: ", round(cg; sigdigits = 3)))
        heatmap!(p, Ds, ϵs, mat'; label = "")
        PltU.sfig(p, 
            plotsdir(simid, "D_ϵ_heatmap", (;cg), ".png")
        )
    end # for cg in cgs
end