# Adapted from:
# SOSDEMO1 --- Sum of Squares Test
# Section 3.1 of SOSTOOLS User's Manual

facts("SOSDEMO1") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x y

  m = Model(solver = solver)

  # Example 2.4 of
  # P. Parrilo and A. Jadbabaie
  # "Approximation of the joint spectral radius using sum of squares."
  # Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402

  p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4

  soscon = @polyconstraint m p >= 0

  status = solve(m)

  @fact status --> :Optimal

  q = getslack(soscon)
  Q = getmat(q)
  @fact issymmetric(Q) --> true
  @fact isapprox(Q[1, 1], 2, rtol=1e-5) --> true
  @fact isapprox(Q[1, 2], 1, rtol=1e-5) --> true
  @fact isapprox(Q[3, 3], 5, rtol=1e-5) --> true
  @fact abs(Q[2, 3]) < 1e-5 --> true
  @fact isapprox(Q[2, 2] + 2Q[1, 3], -1, rtol=1e-5) --> true
  sosdec = SOSDecomposition(q)
  @fact isapprox(sum(sosdec.ps.^2), p; rtol=1e-4) --> true

  M = Model(solver = solver)

  p = 4*x^4*y^6 + x^2 - x*y^2 + y^2

  soscon = @polyconstraint M p >= 0

  status = solve(M)

  @fact status --> :Optimal

  # p should be MatPolynomial([1, 0, -1/2, 0, -1, 1, 0, -2/3, 0, 4/3, 0, 0, 2, 0, 4], [y, x, x*y, x*y^2, x^2*y^3])
end; end; end
