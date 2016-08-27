# Adapted from:
# SOSDEMO6 --- MAX CUT
# Section 3.6 of SOSTOOLS User's Manual

facts("SOSDEMO6") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x1 x2 x3 x4 x5
  x = [x1, x2, x3, x4, x5]

  # Number of cuts
  f = 2.5 - 0.5*x1*x2 - 0.5*x2*x3 - 0.5*x3*x4 - 0.5*x4*x5 - 0.5*x5*x1

  # Boolean constraints
  bc = [x1^2 - 1, x2^2 - 1, x3^2 - 1, x4^2 - 1, x5^2 - 1]

  for (gamma, expected) in [(3.9, :Infeasible), (4, :Optimal)]

    m = JuMP.Model(solver = solver)

    Z = monomials(x, 0:1)
    p = Vector{Any}(6)
    @SOSvariable m p1 >= 0 Z

    Z = monomials(x, 0:2)
    p = Vector{VecPolynomial{Variable}}(5)
    for i in 1:5
      @SOSvariable m tmp Z
      p[i] = tmp
    end

    @SOSconstraint m p1*(gamma-f) + dot(p, bc) >= (gamma-f)^2

    status = solve(m)

    @fact status --> expected
  end
end; end; end
