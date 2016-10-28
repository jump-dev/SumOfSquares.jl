# Adapted from:
# SOSDEMO3 --- Bound on Global Extremum
# Section 3.3 of SOSTOOLS User's Manual

facts("SOSDEMO3") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x1 x2

  m = SOSModel(solver = solver)

  @variable m γ

  # Constraint : r(x)*(f(x) - gam) >= 0
  # f(x) is the Goldstein-Price function
  f1 = x1+x2+1
  f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
  f3 = 2*x1-3*x2
  f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2

  f = (1+f1^2*f2)*(30+f3^2*f4)

  @polyconstraint m f >= γ

  @objective m Max γ

  status = solve(m)

  @fact status --> :Optimal

  @fact getobjectivevalue(m) --> roughly(3; rtol=1e-3)
end; end; end
