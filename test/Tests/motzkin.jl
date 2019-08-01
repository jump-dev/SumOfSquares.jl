using Test
using SumOfSquares
using DynamicPolynomials

function motzkin_test(optimizer, config::MOIT.TestConfig)
    @polyvar x y
    # Motzkin polynomial
    p = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2

    model = _model(optimizer)
    # We want to write `p ≥ 0` instead of `p in SOSCone()` so we need to set
    # the polymodule to `SumOfSquares` to tell `PolyJuMP` to interpret
    # `p ≥ 0` as a Sum-of-Squares constraint.
    PolyJuMP.setpolymodule!(model, SumOfSquares)
    con_ref = @constraint(model, p ≥ 0)

    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.INFEASIBLE
    @test JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE

    q = (x^2 + y^2) * p
    delete(model, con_ref)
    @constraint(model, q >= 0)

    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == config.optimal_status
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
end
sd_tests["motzkin"] = motzkin_test
