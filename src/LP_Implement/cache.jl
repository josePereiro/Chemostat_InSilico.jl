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

    wl = ReentrantLock() # write lock

    # ranges
    vatp_range, vg_ranges = vatpvg_ranges(M; vatp_margin, vg_margin)
    i_vatp_range = vatp_range |> enumerate |> collect
    Nvatpvg = sum(length.(values(vg_ranges)))
    Nvatp = length(vatp_range)

    # Caching
    cache = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    prog = Progress(Nvatp; desc = "Caching (N = $Nvatpvg)  ... ", dt = 0.5)
    c = 0
    @threads for (vatpi, vatp) in i_vatp_range

        lnet = deepcopy(M.net)

        lock(wl) do
            cache[vatp] = Dict{Float64, Vector{Float64}}()
            c += 1
        end

        vg_range = vg_ranges[vatpi]
        for vg in vg_range
            fixxing(lnet, M.vatp_idx, vatp) do 
                fixxing(lnet, M.vg_idx, vg) do 
                    sol = fba(lnet, M.obj_idx)
                    cache[vatp][vg] = sol            
                end
            end
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
