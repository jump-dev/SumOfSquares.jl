# Adapted from:
# SOSDEMO10 --- Set containment
# Section 3.10 of SOSTOOLS User's Manual

function sosdemo10_test(optimizer, config::MOI.Test.Config)
    @polyvar x[1:2]

    ε = 1e-6
    p = x[1]^2 + x[2]^2
    γ = 1
    g0 = 2 * x[1]
    θ = 1

    model = _model(optimizer)
    PolyJuMP.setpolymodule!(model, SumOfSquares)

    # FIXME s should be sos ?
    # in SOSTools doc it is said to be SOS
    # but in the demo it is not constrained so
    Z = monomials(x, 0:4)
    @variable(model, s, Poly(Z))

    Z = monomials(x, 2:3)
    @variable(model, g1, Poly(Z))

    Sc = [
        θ^2-s*(γ-p) g0+g1
        g0+g1 1
    ]

    @constraint(model, Matrix(ε * I, 2, 2) <= Sc, PSDCone())

    JuMP.optimize!(model)

    # Program is feasible, that is, the set
    # { x |((g0 + g1) + θ)(θ - (g0 + g1)) >=0 }
    # contains the set
    # { x | p <= gamma }

    # ALMOST_OPTIMAL and NEARLY_FEASIBLE_POINT for CSDP
    @test JuMP.termination_status(model) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
    @test JuMP.primal_status(model) in
          [MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT]

    #@show JuMP.value(s)
    #@show JuMP.value(g1)
end
sd_tests["sosdemo10"] = sosdemo10_test
