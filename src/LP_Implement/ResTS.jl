struct ResTS
    sg_ts::Vector{Float64}
    sl_ts::Vector{Float64}
    D_ts::Vector{Float64}
    X_ts::Vector{Float64}
end

ResTS() = ResTS([],[],[],[])

function Base.push!(ts::ResTS, M::SimModel)
    push!(ts.sg_ts, M.sg)
    push!(ts.sl_ts, M.sl)
    push!(ts.D_ts, M.D)
    push!(ts.X_ts, M.X)
end
