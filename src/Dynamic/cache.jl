## ---------------------------------------------------------  
_vgvatp_cache_file(M::SimModel, marginf) = joinpath(cachedir(), 
    UJL.mysavename("cache_", "jld"; M.δvatp, M.δvg, marginf)
)

## ---------------------------------------------------------  
# TODO: Runge-Kutta damping
function vgvatp_cache(M::SimModel; marginf::Float64 = 0.1, 
        up_frec = 10)::Dict{Float64, Dict{Float64, Vector{Float64}}}

    # check cache
    cfile = _vgvatp_cache_file(M, marginf)
    isfile(cfile) && return deserialize(cfile)
    
    # setup
    M = deepcopy(M)
    Δb = abs.(M.net.ub .- M.net.lb)
    M.net.ub .+= Δb .* marginf
    M.net.lb .-= Δb .* marginf

    # threading
    nth = nthreads()
    wl = ReentrantLock() # write lock
    net_pool = [deepcopy(M.net) for th in 1:nth]

    # ranges
    vatp_margin = abs(5 * (10.0^(-M.δvatp))) # To get 0 + ϵ values 
    vg_margin = abs(5 * (10.0^(-M.δvg)))
    vatp_range, vg_ranges = vatpvg_ranges(M; vatp_margin, vg_margin)
    i_vatp_range = vatp_range |> enumerate |> collect
    Nvatpvg = sum(length.(values(vg_ranges)))
    Nvatp = length(vatp_range)

    # Caching
    cache = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    prog = Progress(Nvatp; desc = "Caching (N = $Nvatpvg)  ... ", dt = 0.5)
    c = 0
    @threads for (vatpi, vatp) in i_vatp_range

        thid = threadid()
        net = net_pool[thid]

        lock(wl) do
            cache[vatp] = Dict{Float64, Vector{Float64}}()
            c += 1
        end

        vg_range = vg_ranges[vatpi]
        for vg in vg_range
            fixxing(net, M.vatp_idx, vatp) do 
                fixxing(net, M.vg_idx, vg) do 
                    sol = fba(net, M.obj_idx)
                    cache[vatp][vg] = sol            
                end
            end
        end
        up = rem(c, up_frec) == 0
        up && update!(prog, c; showvalues = [
                    (:vatp_range, string(vatp_range, " len: ", length(vatp_range))),
                    (:vg_range, string(vg_range, " len: ", length(vg_range))),
                ]
            )   
    end
    finish!(prog)
    
    serialize(cfile, cache)
    
    return cache
end
