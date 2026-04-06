using Test
using SumOfSquares
using DynamicPolynomials
using JuMP
import Clarabel
import StarAlgebras as SA
import MultivariateBases as MB
import MultivariatePolynomials as MP
using MultivariateMoments: SymMatrix

const nc_optimizer =
    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

@testset "Noncommutative" begin
    @testset "GramMatrix with NC variables" begin
        @ncpolyvar x y
        # NC polynomial algebra_element and polynomial recovery
        p = x * y + y * x + x^2
        a = MB.algebra_element(p)
        @test MP.polynomial(a) == p
        # NC monomials are distinct: xy ≠ yx
        @test x * y != y * x
        @test length(monomials(p)) == 3
    end

    @testset "SOS with NC: (xy + x^2)^2 [BKP16 Example 2.11]" begin
        @ncpolyvar x y
        p = (x * y + x^2)^2
        model = Model(nc_optimizer)
        cref = @constraint(model, p in SOSCone())
        optimize!(model)
        @test termination_status(model) == OPTIMAL

        # Gram matrix and polynomial recovery
        g = gram_matrix(cref)
        @test polynomial(g) ≈ p atol = 1e-5

        # Newton polytope should select [xy, x^2] as basis
        basis_monos = MB.keys_as_monomials(g.basis)
        @test Set(basis_monos) == Set([x * y, x^2])

        # SOS decomposition should have 1 term
        dec = sos_decomposition(cref, 1e-6)
        @test length(dec.ps) == 1
        # The decomposition q should satisfy q' * q ≈ p
        q = polynomial(dec.ps[1])
        @test q' * q ≈ p atol = 1e-5
    end

    @testset "SOS with NC: large Newton polytope [BKP16 Example 2.2]" begin
        @ncpolyvar x y
        n = 5
        p = (x + x^n * y^(2n) * x^n)^2
        model = Model(nc_optimizer)
        cref = @constraint(model, p in SOSCone())
        optimize!(model)
        @test termination_status(model) == OPTIMAL

        # Newton chip should give only 2 basis monomials
        g = gram_matrix(cref)
        basis_monos = MB.keys_as_monomials(g.basis)
        @test length(basis_monos) == 2
        @test Set(basis_monos) == Set([x, x^n * y^(2n) * x^n])
        @test polynomial(g) ≈ p atol = 1e-4
    end

    @testset "Commutator squared is SOS" begin
        @ncpolyvar x y
        comm = x * y - y * x
        p = comm^2  # = xyxy + yxyx - xy^2x - yx^2y
        model = Model(nc_optimizer)
        cref = @constraint(model, p in SOSCone())
        optimize!(model)
        @test termination_status(model) == OPTIMAL
        @test polynomial(gram_matrix(cref)) ≈ p atol = 1e-4
    end

    @testset "Higher-degree commutator squared" begin
        @ncpolyvar x y
        p = (x^2 * y - y * x^2)^2
        model = Model(nc_optimizer)
        cref = @constraint(model, p in SOSCone())
        optimize!(model)
        @test termination_status(model) == OPTIMAL
        dec = sos_decomposition(cref, 1e-6)
        @test length(dec.ps) >= 1
        @test polynomial(gram_matrix(cref)) ≈ p atol = 1e-4
    end

    @testset "Commutator is not SOS" begin
        @ncpolyvar x y
        model = Model(nc_optimizer)
        @constraint(model, x * y - y * x in SOSCone())
        optimize!(model)
        @test termination_status(model) in [INFEASIBLE, INFEASIBLE_OR_UNBOUNDED]
    end

    @testset "NC SOS optimization" begin
        @ncpolyvar x y
        # max t s.t. x^2 + y^2 - t(xy + yx) is SOS
        # Gram matrix Q = [[1, -t], [-t, 1]] on basis [x, y]
        # PSD iff t <= 1
        model = Model(nc_optimizer)
        @variable(model, t)
        @constraint(model, x^2 + y^2 - t * (x * y + y * x) in SOSCone())
        @objective(model, Max, t)
        optimize!(model)
        @test termination_status(model) == OPTIMAL
        @test value(t) ≈ 1.0 atol = 1e-4
    end

    @testset "NC quartic optimization" begin
        @ncpolyvar x y
        # max t s.t. x^4 + y^4 + xyxy + yxyx - t(x^4 + x^2y^2 + y^2x^2 + y^4) SOS
        # This tests optimization with NC degree-4 polynomials
        model = Model(nc_optimizer)
        @variable(model, t)
        p = x^4 + y^4 + x * y * x * y + y * x * y * x
        q = x^4 + x^2 * y^2 + y^2 * x^2 + y^4
        @constraint(model, p - t * q in SOSCone())
        @objective(model, Max, t)
        optimize!(model)
        @test termination_status(model) == OPTIMAL
        @test value(t) ≈ 0.5 atol = 1e-4
    end

    @testset "Hermitian SOS" begin
        @ncpolyvar x y
        p = (x + im * y) * (x - im * y)  # = x^2 + y^2 + i(yx - xy)
        model = Model(nc_optimizer)
        cone =
            NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
        cref = @constraint(model, p in cone)
        optimize!(model)
        @test termination_status(model) == OPTIMAL
        dec = sos_decomposition(cref, 1e-6)
        @test length(dec.ps) == 1
        # Verify decomposition: p* * p should give original
        q = polynomial(dec.ps[1])
        @test q' * q ≈ p atol = 1e-5
    end

    @testset "NC x^2y^2 + y^2x^2 is not SOS" begin
        @ncpolyvar x y
        # In NC, x^2*y^2 + y^2*x^2 is NOT a sum of squares
        p = x^2 * y^2 + y^2 * x^2
        model = Model(nc_optimizer)
        @constraint(model, p in SOSCone())
        optimize!(model)
        @test termination_status(model) in [INFEASIBLE, INFEASIBLE_OR_UNBOUNDED]
    end

    @testset "NC Newton polytope" begin
        @ncpolyvar a b
        uni = Certificate.NewtonDegreeBounds(tuple())
        # Test from certificate.jl
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [a^4, a^3 * b, a * b * a^2, a * b * a * b],
            uni,
        ) == [a * b, a^2]
        # With NewtonFilter
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [
                a^2,
                a^10 * b^20 * a^11,
                a^11 * b^20 * a^10,
                a^10 * b^20 * a^20 * b^20 * a^10,
            ],
            Certificate.NewtonFilter(uni),
        ) == [a, a^10 * b^20 * a^10]
    end
end
