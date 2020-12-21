const DATA_KEY = :DAT
save_data(filename, data) = wsave(filename, Dict(DATA_KEY => data))
load_data(filename) = get(wload(filename), DATA_KEY, nothing)

sci(n) = @sprintf("%0.1e", n)
function mysavename(name, ext = ""; c...)
    d = Dict{Symbol, Any}(c)
    for (k, v) in d
        if v isa AbstractFloat
            d[k] = abs(log10(abs(v))) > 3 ? sci(v) : v
        end
    end
    savename(name, d, ext)
end
