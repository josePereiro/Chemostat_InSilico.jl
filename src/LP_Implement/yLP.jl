## -----------------------------------------------------------------------------------------------
function yLP(S, b, lb, ub, c, d; 
        ϵ = 1e-3, 
        sense = JuMP.MOI.MAX_SENSE,
        solver = GLPK.Optimizer
    )

    lp_model = JuMP.Model(solver)
    M,N = size(S)

    # Variables
    y = JuMP.@variable(lp_model, y[1:N])
    r = JuMP.@variable(lp_model, r)

    # Constraints
    JuMP.@constraint(lp_model, S * y - r .* b .== 0.0)
    JuMP.@constraint(lp_model, d' * y == 1.0)
    JuMP.@constraint(lp_model, y - r .* lb .>= 0.0)
    JuMP.@constraint(lp_model, y - r .* ub .<= 0.0)
    JuMP.@constraint(lp_model, r >= ϵ)

    # objective
    JuMP.@objective(lp_model, sense, c' * y)

    # optimize
    JuMP.optimize!(lp_model)

    yval = JuMP.value.(y)
    rval = JuMP.value(r)
    v = yval ./ rval # sol
    y = (c' * v) / (d' * v) # yield

    status = JuMP.termination_status(lp_model)
    status, v, y
end
yLP(net, d; kwargs...) = yLP(net.S, net.b, net.lb, net.ub, net.c, d; kwargs...)