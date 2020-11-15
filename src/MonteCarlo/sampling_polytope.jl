# Generate N rand cell 
function generate_random_cells(p, n = 1; 
        tries = 5000, verbose = false, 
        threading_th = Int(1e6)
    )

    vglb, vgub = vg_global_min(p), vg_global_max(p)
    vatplb, vatpub = vatp_global_min(p), vatp_global_max(p)    
    
    chuncks = get_chuncks(1:n, nthreads(); th = threading_th)
    rvatps, rvgs = Vector{Float64}(undef, n), Vector{Float64}(undef, n)

    verbose && (prog = Progress(n; desc = "Generating -t$(min(length(chuncks), nthreads()))... ", dt = 0.5); 
        count = zeros(Int, nthreads()))
    @threads for chunck in chuncks
        for i in chunck
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