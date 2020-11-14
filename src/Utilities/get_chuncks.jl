# get a collection and returns an array of chuncks
function get_chuncks(r, n; th = Int(1e4))
    length(r) <= th && return [r]
    chunck_len = length(r) ÷ n + 1
    return [r[(i*chunck_len + 1):min(length(r), (i*chunck_len + chunck_len))] for i in 0:(n - 1)]
end
