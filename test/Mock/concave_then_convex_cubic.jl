config = MOI.Test.Config()
vals = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        vals,
        (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) =>
            zeros(26),
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(2), zeros(2)],
    )
end
for mock in mocks(optimize!)
    Tests.dsos_concave_then_convex_cubic_test(mock, config)
end
function optimize_bridged!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        vals,
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(2), zeros(2)],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [zeros(3) for i in 1:8],
    )
end
function optimize_cached!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        vals,
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(2), zeros(2)],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [zeros(3) for i in 1:8],
    )
end
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sdsos_concave_then_convex_cubic_test(mock, config)
end
function optimize_bridged!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        vals,
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
        (MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle) =>
            [zeros(6), zeros(6)],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(2), zeros(2)],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [zeros(3), zeros(3)],
    )
end
function optimize_cached!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        vals,
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
        (MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle) =>
            [zeros(6), zeros(6)],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(2), zeros(2)],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [zeros(3), zeros(3)],
    )
end
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sos_concave_then_convex_cubic_test(mock, config)
end
