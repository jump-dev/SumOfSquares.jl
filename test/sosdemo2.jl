# Adapted from:
# SOSDEMO2 --- Lyapunov Function Search
# Section 3.2 of SOSTOOLS User's Manual

using Calculus

facts("SOSDEMO2") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x[1:3]

  # Constructing the vector field dx/dt = f
  f = [-x[1]^3-x[1]*x[3]^2,
      -x[2]-x[1]^2*x[2],
      -x[3]+3*x[1]^2*x[3]-3*x[3]/(x[3]^2+1)]

  m = JuMP.Model(solver = solver)

  # The Lyapunov function V(x):
  Z = x.^2
  @SOSvariable m V Z

  @SOSconstraint m V >= sum(x.^2)

  # dV/dx*(x[3]^2+1)*f <= 0
  P = dot(differentiate(V, x), f)*(x[3]^2+1)
  @SOSconstraint m P <= 0

  status = solve(m)

  @fact status --> :Optimal

  @fact removemonomials(getvalue(V), Z) --> zero(VecPolynomial{Float64})
end; end; end
