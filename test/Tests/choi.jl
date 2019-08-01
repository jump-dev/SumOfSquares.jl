using Test
using SumOfSquares
using DynamicPolynomials

"""
    choi_test(optimizer, config::MOIT.TestConfig)

Test using the example of polynomial matrix that is not a Sum-of-Squares matrix
given in [C15] as the following biquadratic form:

F(x; y) = (x_1^2 + 2 x_3^2) * y_1^2
        + (x_2^2 + 2 x_1^2) * y_2^2
        + (x_3^2 + 2 x_2^2) * y_3^2
        - 2 (x_1x_2 y_1y_2 + x_2x_3 y_2y_3 + x_3x_1 y_3y_1)

[C15] Choi, M. D.
*Positive semidefinite biquadratic forms*.
Linear Algebra and its Applications, **1975**, 12(2), 95-100.
"""
function choi_test(optimizer, config::MOIT.TestConfig)
    @polyvar x y z

    model = _model(optimizer)
    PolyJuMP.setpolymodule!(model, SumOfSquares)

    C = [ x^2 + 2y^2 -x*y        -x*z
         -x*y         y^2 + 2z^2 -y*z
         -x*z        -y*z         z^2 + 2x^2]

    @SDconstraint(model, C ⪰ 0)

    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.INFEASIBLE
    @test JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
end
sd_tests["choi_term"] = choi_test
