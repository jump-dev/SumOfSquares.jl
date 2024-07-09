using Test
using SumOfSquares
using DynamicPolynomials

function univariate_sum_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @polyvar x y
    # (x - 1)^2 + (y + 1)^2
    cref = @constraint(
        model,
        x^2 + y^2 + 2(y - x) + 2 in cone,
        sparsity = Sparsity.Variable()
    )

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL

    @test primal_status(model) == MOI.FEASIBLE_POINT

    p = gram_matrix(cref)
    @test p isa SumOfSquares.BlockDiagonalGramMatrix
    @test length(p.blocks) == 2
    @test value_matrix(p.blocks[1]) ≈ [1 -1; -1 1] atol = atol rtol = rtol
    @test p.blocks[1].basis.monomials == [1, x]
    @test value_matrix(p.blocks[2]) ≈ ones(2, 2) atol = atol rtol = rtol
    @test p.blocks[2].basis.monomials == [1, y]

    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace,
        MB.SubBasis{MB.Monomial,monomial_type(x),monomial_vector_type(x)},
        SumOfSquares.Certificate.Sparsity.Ideal{
            Sparsity.Variable,
            SumOfSquares.Certificate.MaxDegree{
                typeof(cone),
                MB.FullBasis{MB.Monomial,monomial_type(x)},
                MB.FullBasis{MB.Monomial,monomial_type(x)},
            },
        },
    }
    @test list_of_constraint_types(model) == [(Vector{AffExpr}, S)]
    return test_delete_bridge(
        model,
        cref,
        0,
        (
            (MOI.VectorOfVariables, MOI.Nonnegatives, 0),
            (MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle, 0),
        ),
    )
end
function sos_univariate_sum_test(optimizer, config)
    return univariate_sum_test(optimizer, config, SOSCone())
end
sd_tests["sos_univariate_sum"] = sos_univariate_sum_test
function sdsos_univariate_sum_test(optimizer, config)
    return univariate_sum_test(optimizer, config, SDSOSCone())
end
soc_tests["sdsos_univariate_sum"] = sdsos_univariate_sum_test
function dsos_univariate_sum_test(optimizer, config)
    return univariate_sum_test(optimizer, config, DSOSCone())
end
linear_tests["dsos_univariate_sum"] = dsos_univariate_sum_test
