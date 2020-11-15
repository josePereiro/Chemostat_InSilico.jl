const DATA_KEY = :DAT
save_data(filename, data) = wsave(filename, Dict(DATA_KEY => data))
load_data(filename) = get(wload(filename), DATA_KEY, nothing)