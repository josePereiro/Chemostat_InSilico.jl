# const PS = DictTree()
function mysavefig(p, pname, dir, id; params...)
    pname = mysavename(pname, "png"; params...)
    fname = joinpath(dir, string(id, "_", pname))
    # PS[pname] = deepcopy(p)
    savefig(p, fname)
    @info "Plotting" fname
end
function mysavefig(ps::Vector, pname, dir, id; 
        layout = _auto_layout(length(ps)), 
        margin = 10, params...
    )
    pname = mysavename(pname, "png"; params...)
    fname = joinpath(dir, string(id, "_", pname))
    grid = make_grid(ps; layout, margin)
    # PS[pname] = deepcopy(grid)
    FileIO.save(fname, grid)
    @info "Plotting" fname
end