# Adapted from:
# SOSDEMO6 --- MAX CUT
# Section 3.6 of SOSTOOLS User's Manual

using Test
using SumOfSquares
using DynamicPolynomials

function maxcut_test(optimizer, config::MOIT.TestConfig)
    atol = config.atol
    rtol = config.rtol

    @polyvar x[1:5]

    # Number of cuts
    f = 2.5 - 0.5*x[1]*x[2] - 0.5*x[2]*x[3] - 0.5*x[3]*x[4] - 0.5*x[4]*x[5] - 0.5*x[5]*x[1]

    # Boolean constraints
    bc = vec(x).^2 .- 1

    model = _model(optimizer)

    Z = monomials(x, 0:1)
    @variable(model, p1, SOSPoly(Z))

    Z = monomials(x, 0:2)
    @variable(model, p[1:5], Poly(Z))

    @testset "with γ=$γ it should be $(feasible ? "feasible" : "infeasible")" for (γ, feasible) in [(3.9, false), (4.1, true)]
        cref = @constraint(model, p1*(γ-f) + dot(p, bc) - (γ-f)^2 in SOSCone())

        JuMP.optimize!(model)

        if feasible
            @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        else
            @test JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        end
        delete(model, cref)
    end
end
sd_tests["maxcut"] = maxcut_test
