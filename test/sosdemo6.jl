# Adapted from:
# SOSDEMO6 --- MAX CUT
# Section 3.6 of SOSTOOLS User's Manual

facts("SOSDEMO6") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x[1:5]

  # Number of cuts
  f = 2.5 - 0.5*x[1]*x[2] - 0.5*x[2]*x[3] - 0.5*x[3]*x[4] - 0.5*x[4]*x[5] - 0.5*x[5]*x[1]

  # Boolean constraints
  bc = x.^2 - 1

  for (gamma, expected) in [(3.9, :Infeasible), (4, :Optimal)]

    m = Model(solver = solver)

    Z = monomials(x, 0:1)
    p = Vector{Any}(6)
    @polyvariable m p1 >= 0 Z

    Z = monomials(x, 0:2)
    p = Vector{VecPolynomial{Variable}}(5)
    for i in 1:5
      @polyvariable m tmp Z
      p[i] = tmp
    end

    @polyconstraint m p1*(gamma-f) + dot(p, bc) >= (gamma-f)^2

    status = solve(m)

    @fact status --> expected
  end
end; end; end
