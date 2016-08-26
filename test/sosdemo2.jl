# Adapted from:
# SOSDEMO2 --- Lyapunov Function Search
# Section 3.2 of SOSTOOLS User's Manual

using Calculus

facts("SOSDEMO2") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x1 x2 x3

  # Constructing the vector field dx/dt = f
  f = [-x1^3-x1*x3^2,
      -x2-x1^2*x2,
      -x3+3*x1^2*x3-3*x3/(x3^2+1)]

  m = JuMP.Model(solver = solver)

  # The Lyapunov function V(x):
  @SOSvariable m V [x1^2, x2^2, x3^2]

  @SOSconstraint m V >= x1^2+x2^2+x3^2

  # dV/dx*(x3^2+1)*f <= 0
  x = [x1, x2, x3]
  P = dot(differentiate(V, x), f)*(x3^2+1)
  @SOSconstraint m P <= 0

  status = solve(m)

  @fact status --> :Optimal

  @fact removemonomials(getvalue(V), [x1^2, x2^2, x3^2]) --> zero(VecPolynomial{Float64})
end; end; end
