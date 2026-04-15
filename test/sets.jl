using Test
using SumOfSquares
import MultivariateBases as MB
import MultivariatePolynomials as MP
import StarAlgebras as SA
import MathOptInterface as MOI
using DynamicPolynomials
using LinearAlgebra

@testset "ScaledMonomial gram matrix" begin
    @polyvar x y

    # For gram basis [y, x] the polynomial is p = b^T Q b.
    # With Q = I (identity), p should be the same regardless of basis.

    # --- Monomial gram basis ---
    # b = [y, x], Q = I ⟹ p = y² + 2·1·xy + x² = y² + 2xy + x²
    # Vectorized Q = [Q₁₁, Q₁₂, Q₂₂] = [1, 1, 1] for the monomial basis
    p_expected = y^2 + 2x * y + x^2

    Q_mono = [1.0, 1.0, 1.0]
    gb_mono = MB.SubBasis{MB.Monomial}([y, x])
    gram_mono = SumOfSquares.build_gram_matrix(
        Q_mono,
        gb_mono,
        MOI.PositiveSemidefiniteConeTriangle,
        Float64,
    )
    p_mono = MP.polynomial(gram_mono)
    @test p_mono ≈ p_expected

    # --- ScaledMonomial gram basis ---
    # With the correct MStruct, ŷ·x̂ has coefficient 1/√2 in the algebra,
    # so b̃^T Q̃ b̃ = Q̃₁₁ ŷ² + 2·(1/√2)·Q̃₁₂ ŝ(xy) + Q̃₂₂ x̂²
    # For Q̃ = I: p = ŷ² + √2·ŝ(xy) + x̂² = y² + √2·√2·xy + x² = y² + 2xy + x²
    Q_scaled = [1.0, 1.0, 1.0]
    gb_scaled = MB.SubBasis{MB.ScaledMonomial}([y, x])
    gram_scaled = SumOfSquares.build_gram_matrix(
        Q_scaled,
        gb_scaled,
        MOI.PositiveSemidefiniteConeTriangle,
        Float64,
    )
    p_scaled = MP.polynomial(gram_scaled)
    @test p_scaled ≈ p_expected
end
