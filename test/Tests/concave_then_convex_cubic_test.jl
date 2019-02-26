using Test
using SumOfSquares
using DynamicPolynomials

function concave_then_convex_cubic_test(optimizer, config::MOIT.TestConfig)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @polyvar x
    @variable(model, p, Poly(monomials(x, 0:3)))
    @constraint(model,  p in SOSConvexCone(), domain = (@set x >= 0))
    @constraint(model, -p in SOSConvexCone(), domain = (@set x <= 0))
    @constraint(model, p(x =>  2) ==  8)
    @constraint(model, p(x => -2) == -8)
    @constraint(model, p(x =>  1) ==  1)
    @constraint(model, p(x => -1) == -1)

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(p) ≈ x^3 atol=atol rtol=rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
end
