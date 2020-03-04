using LinearAlgebra, Test
using SumOfSquares
using DynamicPolynomials

"""
    simple_matrix(optimizer, config::MOIT.TestConfig)

Example 3.77 and 3.79 of
Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
Semidefinite optimization and convex algebraic geometry SIAM 2013
"""
function simple_matrix_test(optimizer, config::MOIT.TestConfig)
    @polyvar x
    P = [x^2-2x+2 x; x x^2]

    # Example 3.77
    mat_model = _model(optimizer)
    PolyJuMP.setpolymodule!(mat_model, SumOfSquares)
    mat_cref = @constraint(mat_model, Symmetric(P) in PSDCone())
    JuMP.optimize!(mat_model)
    @test JuMP.termination_status(mat_model) == MOI.OPTIMAL
    @test JuMP.primal_status(mat_model) == MOI.FEASIBLE_POINT
    @test length(gram_matrix(mat_cref).basis.monomials) == 3

    # Example 3.79
    @polyvar y[1:2]
    model = _model(optimizer)
    PolyJuMP.setpolymodule!(model, SumOfSquares)
    cref = @constraint(model, dot(vec(y), P * vec(y)) >= 0)
    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test length(gram_matrix(cref).basis.monomials) == 3
end
sd_tests["simple_matrix"] = simple_matrix_test
