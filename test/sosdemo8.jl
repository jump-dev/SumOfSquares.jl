# Adapted from:
# SOSDEMO8 --- Bounds in Probability
# Section 3.8 of SOSTOOLS User's Manual

facts("SOSDEMO8") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x

  # The probability adds up to one.
  m0 = 1

  # Mean
  m1  = 1

  # Variance
  sig = 1/2

  # E(x^2)
  m2 = sig^2+m1^2

  # Support of the random variable
  R = [0,5]

  # Event whose probability we want to bound
  E = [4,5]

  m = Model()

  @variable m a
  @variable m b
  @variable m c

  P = a + b*x + c*x^2

  # Nonnegative on the support
  @polyconstrant m P >= 0 R

  # Greater than one on the event
  @polyconstrant m P >= 1 E

  # The bound
  bnd =  a * m0 + b * m1 + c * m2

  @objective m Min bnd

  status = solve(m)

  BND = getvalue(bnd,16)
  PP = getvalue(P)
end; end; end
