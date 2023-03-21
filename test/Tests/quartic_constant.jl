using Test
using SumOfSquares
using DynamicPolynomials

function quartic_constant_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
)
    atol = config.atol
    rtol = config.rtol

    @polyvar x
    p = x^4 + 4

    model = _model(optimizer)
    @variable(model, γ)

    cref =
        @constraint(model, p - γ in cone, basis = FixedPolynomialBasis([x^2]))

    @objective(model, Max, γ)

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL

    @test primal_status(model) == MOI.FEASIBLE_POINT

    @test objective_value(model) ≈ 4 atol = atol rtol = rtol

    p = gram_matrix(cref)
    @test p isa SumOfSquares.GramMatrix
    @test getmat(p) ≈ ones(1, 1) atol = atol rtol = rtol
    @test p.basis isa FixedPolynomialBasis
    @test p.basis.polynomials == [x^2]

    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace,
        Monomial{true},
        MonomialVector{true},
        SumOfSquares.Certificate.FixedBasis{typeof(cone),typeof(p.basis)},
    }
    @test list_of_constraint_types(model) == [(Vector{AffExpr}, S)]
    return test_delete_bridge(
        model,
        cref,
        1,
        (
            (MOI.VectorOfVariables, MOI.Nonnegatives, 0),
            (MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle, 0),
        ),
    )
end
function sos_quartic_constant_test(optimizer, config)
    return quartic_constant_test(optimizer, config, SOSCone())
end
sd_tests["sos_quartic_constant"] = sos_quartic_constant_test
function sdsos_quartic_constant_test(optimizer, config)
    return quartic_constant_test(optimizer, config, SDSOSCone())
end
soc_tests["sdsos_quartic_constant"] = sdsos_quartic_constant_test
function dsos_quartic_constant_test(optimizer, config)
    return quartic_constant_test(optimizer, config, DSOSCone())
end
linear_tests["dsos_quartic_constant"] = dsos_quartic_constant_test
