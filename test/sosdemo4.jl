# Adapted from:
# SOSDEMO4 --- Matrix Copositivity
# Section 3.4 of SOSTOOLS User's Manual

facts("SOSDEMO4") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x1 x2 x3 x4 x5
  x = [x1, x2, x3, x4, x5]

  # The matrix under consideration
  J = [1 -1  1  1 -1;
      -1  1 -1  1  1;
       1 -1  1 -1  1;
       1  1 -1  1 -1;
      -1  1  1 -1  1]

  xs = x.^2
  xsJxs = dot(xs, J*xs)
  r = x1^2 + x2^2 + x3^2 + x4^2 + x5^2

  m0 = JuMP.Model(solver = solver)
  @SOSconstraint m0 xsJxs >= 0
  status = solve(m0)
  @fact status --> :Infeasible

  m1 = JuMP.Model(solver = solver)
  @SOSconstraint m1 r*xsJxs >= 0
  status = solve(m1)
  @fact status --> :Optimal
  # Program is feasible. The matrix J is copositive

end; end; end
