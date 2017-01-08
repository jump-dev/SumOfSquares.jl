# Adapted from:
# SOSDEMO7 --- Chebyshev polynomials
# Section 3.7 of SOSTOOLS User's Manual

facts("SOSDEMO7") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    if !isscs(solver)
        ndeg = 8   # Degree of Chebyshev polynomial

        @polyvar x

        Z = monomials([x], 0:ndeg-1)

        m = Model(solver = solver)

        @variable m γ
        @polyvariable m p1 Z

        p = p1 + γ * x^ndeg # the leading coeff of p is γ

        @polyconstraint(m, p <= 1, domain = x >= -1 && x <= 1)
        @polyconstraint(m, p >= -1, domain = x >= -1 && x <= 1)

        @objective m Max γ

        status = solve(m)

        @fact status --> :Optimal

        @fact isapprox(getvalue(p), 128x^8 - 256x^6 + 160x^4 - 32x^2 + 1) --> true
        @fact isapprox(getvalue(γ), 128) --> true
    end
end; end; end
