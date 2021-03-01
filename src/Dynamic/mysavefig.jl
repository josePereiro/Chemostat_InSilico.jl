# const PS = DictTree()
function mysavefig(p, pname, dir, id; params...)
    pname = UJL.mysavename(pname, "png"; params...)
    fname = joinpath(dir, string(id, "_", pname))
    savefig(p, fname)
    @info "Plotting" fname
end
function mysavefig(ps::Vector, pname, dir, id; 
        layout = UJL._auto_layout(length(ps)), 
        margin = 10, params...
    )
    pname = UJL.mysavename(pname, "png"; params...)
    fname = joinpath(dir, string(id, "_", pname))
    grid = UJL.make_grid(ps; layout, margin)
    FileIO.save(fname, grid)
    @info "Plotting" fname
end