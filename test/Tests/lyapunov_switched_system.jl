# Computes a polynomial Lyapunov for the switched system from Example 2.8 of
# [PJ08].
#   x_{k+1} = A_σ x_k, σ ∈ {1, 2}
# Where
#   A1 = [1 0   A2 = [0  1
#         1 0]        0 -1]
# There is a quartic Lyapunov function for the system but only a quadratic for
# the system with A1' = A1/√2, A2' = A2/√2.
# See https://github.com/blegat/SwitchOnSafety.jl for computing Lyapunov
# functions for switched systems.
#
# [PJ08] P. Parrilo and A. Jadbabaie
# Approximation of the joint spectral radius using sum of squares
# Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402
# Inspired from the construction of:
# Ando, T. and Shih, M.-h.
# Simultaneous Contractibility.
# SIAM Journal on Matrix Analysis & Applications, 1998, 19, 487

using LinearAlgebra, Test
using SumOfSquares
using DynamicPolynomials

# Search for polynomial Lyapunov functions of degree `2d` for `A1/γ` and `A2/γ`
function lyapunov_switched_system_test(
    optimizer, config::MOIT.TestConfig, degree::Int, γ, feasible::Bool)

    atol = config.atol
    rtol = config.rtol

    A1 = [1 0
          1 0]
    A2 = [0  1
          0 -1]

    model = _model(optimizer)

    @polyvar x[1:2]

    # p is should be strictly positive so we create a nonnegative `p`
    # and sum it with a strictly positive `q`.
    # Since the problem is homogeneous (i.e. given any `λ > 0`, `p` is
    # feasible iff `λp` is feasible), this is wlog.
    p0 = @variable(model, variable_type=SOSPoly(monomials(x, degree)))
    q = GramMatrix(SOSDecomposition(x.^degree))

    # Keep `p` in a `GramMatrix` form while `q + p0` would transform it to
    # a polynomial. It is not mandatory to keep it in its `GramMatrix` form
    # but it will allow us to test `JuMP.value(::GramMatrix{JuMP.AffExpr})`.
    p = gram_operate(+, q, p0)

    c1 = @constraint(model, p - p(x => (A1 / γ) * vec(x)) in SOSCone())
    c2 = @constraint(model, p - p(x => (A2 / γ) * vec(x)) in SOSCone())

    JuMP.optimize!(model)

    if feasible
        @test JuMP.termination_status(model) == MOI.OPTIMAL
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test all(eigvals(Matrix(getmat(JuMP.value(p0)))) .≥ -atol)
        @test getmat(JuMP.value(p)) ≈ getmat(gram_operate(+, q, JuMP.value(p0))) atol=atol rtol=rtol
    else
        @test JuMP.termination_status(model) == MOI.INFEASIBLE
        @test JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        for (μ1, μ2) in [(JuMP.dual(c1), JuMP.dual(c2)),
                         (moments(c1), moments(c2))]
            # The dual constraint should work on any polynomial.
            # Let's test it with q
            lhs = dot(μ1, q(x => A1 * vec(x))) + dot(μ2, q(x => A2 * vec(x)))
            rhs = dot(μ1, q) + dot(μ2, q)
            @test atol + rtol * max(abs(lhs), abs(rhs)) + lhs >= rhs
        end
    end
end

quadratic_infeasible_lyapunov_switched_system_test(optimizer, config, ε = 2e-1) = lyapunov_switched_system_test(optimizer, config, 1, √2 - ε, false)
sd_tests["quadratic_infeasible_lyapunov_switched_system"] = quadratic_infeasible_lyapunov_switched_system_test
quadratic_feasible_lyapunov_switched_system_test(optimizer,   config, ε = 2e-1) = lyapunov_switched_system_test(optimizer, config, 1, √2 + ε, true)
sd_tests["quadratic_feasible_lyapunov_switched_system"]   = quadratic_feasible_lyapunov_switched_system_test
quartic_infeasible_lyapunov_switched_system_test(optimizer,   config, ε = 2e-1) = lyapunov_switched_system_test(optimizer, config, 2,  1 - ε, false)
sd_tests["quartic_infeasible_lyapunov_switched_system"]   = quartic_infeasible_lyapunov_switched_system_test
quartic_feasible_lyapunov_switched_system_test(optimizer,     config, ε = 2e-1) = lyapunov_switched_system_test(optimizer, config, 2,  1 + ε, true)
sd_tests["quartic_feasible_lyapunov_switched_system"]     = quartic_feasible_lyapunov_switched_system_test
