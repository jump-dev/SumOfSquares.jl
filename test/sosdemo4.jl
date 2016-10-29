# Adapted from:
# SOSDEMO4 --- Matrix Copositivity
# Section 3.4 of SOSTOOLS User's Manual

facts("SOSDEMO4") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x[1:5]

  # The matrix under consideration
  J = [1 -1  1  1 -1;
      -1  1 -1  1  1;
       1 -1  1 -1  1;
       1  1 -1  1 -1;
      -1  1  1 -1  1]

  xs = x.^2
  xsJxs = dot(xs, J*xs)
  r = sum(xs)

  m0 = Model(solver = solver)
  @polyconstraint m0 xsJxs >= 0
  status = solve(m0)
  @fact status --> :Infeasible

  m1 = Model(solver = solver)
  @polyconstraint m1 r*xsJxs >= 0
  status = solve(m1)
  @fact status --> :Optimal
  # Program is feasible. The matrix J is copositive

end; end; end
