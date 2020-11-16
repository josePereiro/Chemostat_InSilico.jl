## ------------------------------------------------------------------
function runMC(; p::Polytope, pinking_fun::Function, 
        ncells::Int, niters::Int, mutr::Float64,
        verbose::Bool = true, threading_th::Int = Int(1e6), 
        feedback_fun::Function = (count, cells_pool) -> false, 
        feedback_frec::Int = max(1, threading_th ÷ 10),
        at_mutate::Function, at_notmutate::Function
    )

    # lock
    # l = ReentrantLock();
    l = SpinLock()

    # Generating pool
    cells_pool = generate_random_cells(p, ncells; verbose, threading_th)
    
    # Chuncks
    chuncks = get_chuncks(1:niters, nthreads(); th = threading_th)

    # Simulating time
    verbose && (prog = Progress(niters; desc = "Simulating -t$(min(length(chuncks), nthreads()))... ", dt = 0.5); 
        count = zeros(Int, nthreads()))
    @threads for chunck in chuncks
        chunck_len = length(chunck)
        for (i, it) in enumerate(chunck)
            pcell = first(pick_cells(pinking_fun, 1, cells_pool))
            ridx = rand(1:ncells)
            cells_pool[ridx] = rand() < mutr ? at_mutate(pcell) : at_notmutate(pcell)
            
            # Feed back
            if i == 1 || i % feedback_frec == 0 || i == chunck_len
                lock(l) do
                    count[threadid()] += feedback_frec
                    ccount = sum(count)
                    feedback_fun(ccount, cells_pool) && return cells_pool
                    verbose && update!(prog, ccount)
                end
            end
        end
    end
    verbose && finish!(prog)
    cells_pool
end