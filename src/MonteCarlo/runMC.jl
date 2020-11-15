## ------------------------------------------------------------------
function generate_dougther(c::Cell)
    vglb, vgub = vg_local_min(c.vatp, c.p), vg_local_max(c.vatp, c.p)
    rvg = rand()*(vgub - vglb) + vglb
    Cell(c.p, c.vatp, rvg)
end

## ------------------------------------------------------------------
function runMC(fedback_fun; p::Polytope, pinking_fun::Function, 
        ncells::Int, niters::Int, mutr::Float64,
        verbose = true, threading_th::Int = Int(1e6)
    )
        
    # Params and tools
    itmutates() = rand() < mutr

    # Generating pool
    cells_pool = generate_random_cells(p, ncells; verbose, threading_th)
    
    # Chuncks
    chuncks = get_chuncks(1:niters, nthreads(); th = threading_th)

    # Simulating time
    verbose && (prog = Progress(niters; desc = "Simulating -t$(min(length(chuncks), nthreads()))... ", dt = 0.5); 
        count = zeros(Int, nthreads()))
    @threads for chunck in chuncks

        for it in chunck
            pcell = first(pick_cells(pinking_fun, 1, cells_pool))
            ridx = rand(1:ncells)
            cells_pool[ridx] = itmutates() ? first(generate_random_cells(p, 1; threading_th = 0)) : generate_dougther(pcell)
            
            # Feed back
            count[threadid()] += 1
            ccount = sum(count)
            fedback_fun(ccount, cells_pool) && return cells_pool

            verbose && update!(prog, ccount)
        end
    end
    verbose && finish!(prog)
    cells_pool
end

runMC(; p::Polytope, pinking_fun::Function, ncells::Int, niters::Int, mutr::Float64, 
        verbose = true, threading_th::Int = Int(1e6)) = 
    runMC((count, cells_pool) -> false; p, pinking_fun, ncells, niters, mutr, verbose, threading_th)