@polyvar x y

m = JuMP.Model(solver = SCSSolver(verbose=false))

p = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2

@SOSconstraint m p >= 0

status = solve(m)

@test status == :Infeasible

M = JuMP.Model(solver = SCSSolver(verbose=false))

q = (x^2 + y^2) * p

@SOSconstraint M q >= 0

status = solve(M)

@test status == :Optimal
