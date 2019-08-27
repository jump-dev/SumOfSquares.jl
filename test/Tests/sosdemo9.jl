# Adapted from:
# SOSDEMO9 --- Matrix SOS decomposition
# Section 3.9 of SOSTOOLS User's Manual

function sosdemo9_test(optimizer, config::MOIT.TestConfig)
    @polyvar x1 x2 x3

    P = [x1^4+x1^2*x2^2+x1^2*x3^2 x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2);
         x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2) x1^2*x2^2+x2^2*x3^2+(x2^2+2*x3^2)^2]

    # Test if P(x1,x2,x3) is an SOS matrix
    model = _model(optimizer)
    PolyJuMP.setpolymodule!(model, SumOfSquares)
    # TODO return H so that P = H.'*H
    @SDconstraint(model, P ⪰ 0)

    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
end
sd_tests["sosdemo9"] = sosdemo9_test
