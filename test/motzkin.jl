facts("Motzkin") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x y

  m = SOSModel(solver = solver)

  p = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2

  @polyconstraint m p >= 0

  status = solve(m)

  @fact status --> :Infeasible

  M = SOSModel(solver = solver)

  q = (x^2 + y^2) * p

  @polyconstraint M q >= 0

  status = solve(M)

  @fact status --> :Optimal
end; end; end
