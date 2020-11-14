function is_inpolytope(vatp, vg, p)
    vglb, vgub = vg_local_min(vatp, p), vg_local_max(vatp, p)
    vatplb, vatpub = vatp_local_min(vglb, p), vatp_local_max(vgub, p)
    vglb <= vg <= vgub && vatplb <= vatp <= vatpub
end

# Generate N rand cell 
function generate_random_cells(p, n = 1; tries = 5000, verbose = false)
    vglb, vgub = vg_global_min(p), vg_global_max(p)
    vatplb, vatpub = vatp_global_min(p), vatp_global_max(p)
    
    verbose && (prog = Progress(n, "Generating -t$(nthreads())... "); count = zeros(Int, nthreads()))
    rvatps, rvgs = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
    @threads for i in collect(1:n)
        found = false
        for t in 1:tries
            rvatps[i] = rand()*(vatpub - vatplb) + vatplb
            rvgs[i] = rand()*(vgub - vglb) + vglb
            is_inpolytope(rvatps[i], rvgs[i], p) &&  (found = true; break)
        end
        if !found
            rvatps[i] = rand()*(vgub - vglb) + vglb
            vatplb, vatpub = vatp_local_min(rvatps[i], p), vatp_local_max(rvatps[i], p)
            rvgs[i] = rand()*(vatpub - vatplb) + vatplb
        end
        verbose && (count[threadid()] += 1; update!(prog, sum(count)))
    end
    verbose && (finish!(prog))
    Cell.([p], rvatps, rvgs)
end


function pick_cells(p::Function, n, cells_pool; 
        verbose = false,
        tries = 500, 
    )
    cells = Vector{Cell}(undef, n)
    for i in 1:n
        cells[i] = rand(cells_pool)
        for t in 1:tries
            rcell = rand(cells_pool)
            p(rcell) && (cells[i] = rcell; break)
        end
    end
    return cells
end