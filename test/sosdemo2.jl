# Adapted from:
# SOSDEMO2 --- Lyapunov Function Search
# Section 3.2 of SOSTOOLS User's Manual

@testset "SOSDEMO2 with $solver" for solver in sdp_solvers
  @polyvar x[1:3]

  # Constructing the vector field dx/dt = f
  f = [-x[1]^3-x[1]*x[3]^2,
      -x[2]-x[1]^2*x[2],
      -x[3]+3*x[1]^2*x[3]-3*x[3]/(x[3]^2+1)]

  m = SOSModel(solver = solver)

  # The Lyapunov function V(x):
  Z = x.^2
  @polyvariable m V Z

  @polyconstraint m V >= sum(x.^2)

  # dV/dx*(x[3]^2+1)*f <= 0
  P = dot(differentiate(V, x), f)*(x[3]^2+1)
  @polyconstraint m P <= 0

  status = solve(m)

  @test status == :Optimal

  @test removemonomials(getvalue(V), Z) == zero(Polynomial{true, Float64})
end
