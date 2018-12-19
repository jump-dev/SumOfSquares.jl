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

@testset "[PJ08] Example 2.8 with $(factory.constructor)" for factory in sdp_factories
    isscs(factory) && continue
    @polyvar x[1:2]
    A1 = [1 0; 1 0]
    A2 = [0 1; 0 -1]
    expected_ub = [√2, 1]
    function testlyap(d, γ, feasible::Bool)
        m = SOSModel(factory)
        @variable m p Poly(monomials(x, 2d))
        # p strictly positive
        q = sum(x.^(2*d))
        @constraint m p >= q
        c1 = @constraint m p(x => A1 * vec(x)) <= γ^(2*d) * p
        c2 = @constraint m p(x => A2 * vec(x)) <= γ^(2*d) * p

        JuMP.optimize!(m)

        if feasible
            @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        else
            @test JuMP.dual_status(m) == MOI.INFEASIBILITY_CERTIFICATE
            μ1 = JuMP.dual(c1)
            μ2 = JuMP.dual(c2)

            # The dual constraint should work on any polynomial.
            # Let's test it with q
            lhs = dot(μ1, q(x => A1 * vec(x))) + dot(μ2, q(x => A2 * vec(x)))
            rhs = dot(μ1, q) + dot(μ2, q)
            @test 1e-6 * max(abs(lhs), abs(rhs)) + lhs >= rhs
        end
    end
    testlyap(1, √2 - 1e-1, false)
    testlyap(1, √2 + 1e-1, true)
    testlyap(2, 1 - 1e-1, false)
    testlyap(2, 1 + 1e-1, true)
end
