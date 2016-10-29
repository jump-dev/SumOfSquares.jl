# Adapted from:
# SOSDEMO7 --- Chebyshev polynomials
# Section 3.7 of SOSTOOLS User's Manual

facts("SOSDEMO7") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  ndeg = 8   # Degree of Chebyshev polynomial

  @polyvar x

  Z = monomials(x, 0:ndeg-1)

  m = Model(solver = solver)

  @polyvariable m p1 Z

  p = p1 + gam * x^ndeg # the leading coeff of p is gam

  @polyconstraint m p <= 1 domain=[x>=-1, x<=1]
  @polyconstraint m p >= -1 domain=[x>=-1, x<=1]

  @objective m Max gam

  status = solve(m)

  @fact status --> :Optimal

  @show getvalue(p)
  @show getvalue(gam)
end; end; end
