function polytope_area(p; 
        vglb = vg_global_min(p), vgub = vg_global_max(p), 
        n::Int = Int(1e9),
        batch_size::Int = Int(1e7), verbose = false)

    batch_size = min(batch_size, n)
    r = range(vglb, vgub; length = n)
    Δvg = step(r)
    vgs = Iterators.Stateful(r)

    area = zeros(nthreads())
    verbose && (prog = Progress(n, "Computing Area -t$(nthreads())... "))
    while length(vgs) > 0
        batch = collect(Iterators.take(vgs, batch_size))
        @threads for vg in batch
            area[threadid()] += Δvatp(vg, p) * Δvg
        end
        verbose && update!(prog, n - length(vgs); showvalues = 
            [
                ("n    ", n),
                ("left ", length(vgs)),
                ("area ", sum(area)),
            ]
        )
    end
    verbose && finish!(prog)
    sum(area)
end
