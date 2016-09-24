# Adapted from:
# SOSDEMO10 --- Set containment
# Section 3.10 of SOSTOOLS User's Manual

facts("SOSDEMO10") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x1 x2
  x = [x1, x2]

  eps = 1e-6
  p = x1^2+x2^2
  gamma = 1
  g0 = 2*x1
  theta = 1

  m = JuMP.Model(solver = solver)

  # FIXME s should be sos ?
  # in SOSTools doc it is said to be SOS
  # but in the demo it is not constrained so
  Z = monomials(x, 0:4)
  @SOSvariable m s Z

  Z = monomials(x, 2:3)
  @SOSvariable m g1 Z

  Sc = [theta^2-s*(gamma-p) g0+g1; g0+g1 1]

  @SOSconstraint m Sc >= eps * eye(2)

  status = solve(m)

  # Program is feasible, { x |((g0+g1) + theta)(theta - (g0+g1)) >=0 } contains { x | p <= gamma }
  @fact status --> :Optimal

  #@show getvalue(s)
  #@show getvalue(g1)
end; end; end
