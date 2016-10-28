# Adapted from:
# SOSDEMO9 --- Matrix SOS decomposition
# Section 3.9 of SOSTOOLS User's Manual

facts("SOSDEMO9") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do

  @polyvar x1 x2 x3

  P = [x1^4+x1^2*x2^2+x1^2*x3^2 x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2);
     x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2) x1^2*x2^2+x2^2*x3^2+(x2^2+2*x3^2)^2]

  # Test if P(x1,x2,x3) is an SOS matrix
  m = SOSModel(solver = solver)
  # TODO return H so that P = H.'*H
  @polyconstraint m P ⪰ 0

  status = solve(m)
  @fact status --> :Optimal
end; end; end
