facts("Monomial selection for certificate") do
  @polyvar x y
  @fact_throws ErrorException getmonomialsforcertificate([x*y, y^2], :Sparse)
  @fact getmonomialsforcertificate([x*y, y^2]) --> MonomialVector([y])
end

facts("Random SOS should be SOS") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x y
  x = [Monomial(), x, y, x^2, y^2, x*y]
  @fact typeof(x) --> Vector{Monomial}
  for i in 1:10
    p = randsos(x)

    m = SOSModel(solver = solver)

    @polyconstraint m p >= 0

    status = solve(m)

    @fact status --> :Optimal
  end
end; end; end
