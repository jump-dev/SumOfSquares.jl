config = MOI.Test.Config(atol=1e-6)
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [2.0, 1.0, 1.0, √2],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [[1.0, 1.0, -√2]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0, -1.0, 1.0]],
    )
end
for mock in mocks(optimize!)
    Tests.sos_univariate_quadratic_test(mock, config)
    Tests.sos_bivariate_quadratic_test(mock, config)
    Tests.sdsos_univariate_quadratic_test(mock, config)
    Tests.sdsos_bivariate_quadratic_test(mock, config)
    Tests.sos_scaled_univariate_quadratic_test(mock, config)
    Tests.sdsos_scaled_univariate_quadratic_test(mock, config)
    Tests.sos_cheby_univariate_quadratic_test(mock, config)
    Tests.sdsos_cheby_univariate_quadratic_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [2.0, 1.0, 1.0, √2],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [[1.0, 1.0, -√2]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0, -√2, 1.0]],
    )
end
for mock in mocks(optimize!)
    Tests.sos_scaled_bivariate_quadratic_test(mock, config)
    Tests.sdsos_scaled_bivariate_quadratic_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [2.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0],
        (
            MOI.VectorOfVariables,
            MOI.PositiveSemidefiniteConeTriangle,
        ) => [[0.0, 0.0, 1.0, 0.0, -1.0, 1.0]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) =>
            [[0.0, 0.0, 0.0, 2.0, -1.0, 2.0]],
    )
end
for mock in mocks(optimize!)
    Tests.sos_cheby_bivariate_quadratic_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        (MOI.FEASIBLE_POINT, [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, √2]),
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, -√2]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) =>
            [[0.0, 0.0, 0.0, 2.0, -1.0, 2.0]],
    )
end
for mock in mocks(optimize!)
    Tests.sdsos_cheby_bivariate_quadratic_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [2.0, 1.0, 1.0, 1.0, 1.0],
        (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) =>
            [0.0, 2.0, 2.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[2.0, -1.0, 0.0]],
    )
end
for mock in mocks(optimize!)
    Tests.dsos_univariate_quadratic_test(mock, config)
    Tests.dsos_bivariate_quadratic_test(mock, config)
    Tests.dsos_scaled_univariate_quadratic_test(mock, config)
    Tests.dsos_cheby_univariate_quadratic_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(mock,
(MOI.FEASIBLE_POINT, [2.0, 1.0, 1.0, 1.0, 1.0]),
(MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) => [0.0, 2.0, 1.0, 1.0],
(MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0, -√2, 1.0]],
)
end
for mock in mocks(optimize!)
    Tests.dsos_scaled_bivariate_quadratic_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(mock,
        (MOI.FEASIBLE_POINT, [2.0, 0.0, 0.0, 1.0, -0.0, 1.0, 1.0, 0.0, 0.0, 1.0]),
        (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) => [0.7624324421943409, 0.7624324421943427, 0.7624324421943439, 0.7624324421943471, 0.0, 2.0, 0.5248650317582692, 1.0, 1.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.5248650317582707, 0.0, 0.0, 1.4751346735031994, -1.0, 1.4751346735032143]],
    )
end
for mock in mocks(optimize!)
    Tests.dsos_cheby_bivariate_quadratic_test(mock, config)
end
