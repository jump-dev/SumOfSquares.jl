# See https://github.com/blegat/SwitchedSystems.jl

# Example 2.8 of
# P. Parrilo and A. Jadbabaie
# Approximation of the joint spectral radius using sum of squares
# Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402
# Inspired from the construction of:
# Ando, T. and Shih, M.-h.
# Simultaneous Contractibility.
# SIAM Journal on Matrix Analysis & Applications, 1998, 19, 487
# The JSR is √2

@testset "[PJ08] Example 2.8 with $solver" for solver in sdp_solvers
    @polyvar x[1:2]
    A1 = [1 0; 1 0]
    A2 = [0 1; 0 -1]
    expected_ub = [√2, 1]
    function testlyap(d, γ, expected_status)
        m = SOSModel(solver = solver)
        X = monomials(x, d)
        @variable m p Poly{false, :Gram}(X)
        # p strictly positive
        q = sum(x.^(2*d))
        @constraint m p >= q
        c1 = @constraint m p(A1 * x, x) <= γ^(2*d) * p
        c2 = @constraint m p(A2 * x, x) <= γ^(2*d) * p

        @test expected_status == solve(m; suppress_warnings=true)

        if expected_status == :Infeasible
            μ1 = getdual(c1)
            μ2 = getdual(c2)

            # The dual constraint should work on any polynomial.
            # Let's test it with q
            lhs = dot(μ1, q(A1 * x, x)) + dot(μ2, q(A2 * x, x))
            rhs = dot(μ1, q) + dot(μ2, q)
            @test 1e-6 * max(abs(lhs), abs(rhs)) + lhs >= rhs
        end
    end
    testlyap(1, √2 - 1e-4, :Infeasible)
    testlyap(1, √2 + 1e-3, :Optimal)
    testlyap(2, 1 - 1e-3, :Infeasible)
    testlyap(2, 1 + 1e-2, :Optimal)
end
