## ---------------------------------------------------------  
cache_file(M::SimModel, marginf) = joinpath(CACHE_DIR, 
    mysavename("cache_", "jld"; M.θvatp, M.θvg, marginf)
)

## ---------------------------------------------------------  
# TODO: Runge-Kutta damping
function vgvatp_cache(M::SimModel; marginf::Int = 50)::Dict{Float64, Dict{Float64, Vector{Float64}}}

    vatp_margin = abs(marginf * (10.0^(-M.θvatp)))
    vg_margin = abs(marginf * (10.0^(-M.θvg)))

    # check cache
    cfile = cache_file(M, marginf)
    isfile(cfile) && return deserialize(cfile)

    # Open network
    M = deepcopy(M)

    # ranges
    vatp_range, vg_ranges = vatpvg_ranges(M; vatp_margin, vg_margin)
    N = sum(length.(values(vg_ranges)))

    # Caching
    cache = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    prog = Progress(N; desc = "Caching (N = $N)  ... ", dt = 0.5)
    c = 0
    for (vatpi, vatp) in enumerate(vatp_range)
        lnet = deepcopy(M.net)
        cache[vatp] = Dict{Float64, Vector{Float64}}()
        vg_range = vg_ranges[vatpi]
        for vg in vg_range
            fixxing(lnet, M.vatp_idx, vatp) do 
                fixxing(lnet, M.vg_idx, vg) do 
                    sol = fba(lnet, M.obj_idx)
                    cache[vatp][vg] = sol            
                end
            end
            c += 1
        end
        update!(prog, c; showvalues = [
                    (:vatp_range, string(vatp_range, " len: ", length(vatp_range))),
                    (:vg_range, string(vg_range, " len: ", length(vg_range))),
                ]
            )   
    end
    finish!(prog)
    
    serialize(cfile, cache)
    
    return cache
end
