config = MOI.Test.TestConfig()
optimize!(mock) = MOIU.mock_optimize!(mock, [2.0, 1.0, 1.0, √2],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[1.0, 1.0, -√2]],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0, -1.0, 1.0]])
for mock in mocks(optimize!)
    Tests.sos_univariate_quadratic_test(mock, config)
    Tests.sos_bivariate_quadratic_test(mock, config)
    Tests.sdsos_univariate_quadratic_test(mock, config)
    Tests.sdsos_bivariate_quadratic_test(mock, config)
end
optimize!(mock) = MOIU.mock_optimize!(mock, [2.0, 1.0, 1.0, 1.0, 1.0],
    (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) => [0.0, 2.0, 2.0, 0.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[2.0, -1.0, 0.0]])
for mock in mocks(optimize!)
    Tests.dsos_univariate_quadratic_test(mock, config)
    Tests.dsos_bivariate_quadratic_test(mock, config)
end
