function generate_random_point(p::Polytope, 
        vatpGL = vatpL(p), ΔvatpG = Δvatp(p),
        vgGL = vgL(p), ΔvgG = Δvg(p); 
        tries = 500
    )
    # try with given coordinates
    for t in 1:tries
        rvatp = rand() * ΔvatpG + vatpGL
        rvg = rand() * ΔvgG + vgGL
        is_inpolytope(rvatp, rvg, p) && return rvatp, rvg
    end

    # find for sure
    rvatp = rand_vatp(p)
    rvg = rand_vg(rvatp, p)

    return rvatp, rvg
end

# Generate N rand cell 
function generate_random_cells(p::Polytope, n::Int = 1; 
        tries::Int = 500, verbose::Bool = false, 
        threading_th::Int = Int(1e6)
    )

    # globals
    vatpGL, ΔvatpG = vatpL(p), Δvatp(p)
    vgGL, ΔvgG = vgL(p), Δvg(p)

    if n == 1
        rvatp, rvg = generate_random_point(p, vatpGL, ΔvatpG, vgGL, ΔvgG)
        return Cell[Cell(rvatp, rvg, p)]
    end
    
    chuncks = get_chuncks(1:n, nthreads(); th = threading_th)
    rvatps, rvgs = Vector{Float64}(undef, n), Vector{Float64}(undef, n)

    verbose && (prog = Progress(n; desc = "Generating -t$(min(length(chuncks), nthreads()))... ", dt = 0.5); 
        count = zeros(Int, nthreads()))
    @threads for chunck in chuncks
        for i in chunck
            rvatps[i], rvgs[i] = generate_random_point(p, vatpGL, ΔvatpG, vgGL, ΔvgG)
            verbose && (count[threadid()] += 1; update!(prog, sum(count)))
        end
    end
    verbose && (finish!(prog))
    Cell.([p], rvatps, rvgs)
end


function pick_cells(p::Function, n::Int, cells_pool::Vector{Cell}; 
        verbose::Bool = false,
        tries::Int = 500, 
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

function generate_similar_cell(c::Cell, dvatp)
    svatpL = max(c.vatp - dvatp, vatpL(c.p))
    svatpU = min(c.vatp + dvatp, vatpU(c.p))
    rvatp = rand() * (svatpU - svatpL) + svatpL
    rvg = rand_vg(rvatp, c.p)
    Cell(c.p, rvatp, rvg)
end
