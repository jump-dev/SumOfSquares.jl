config = MOI.Test.Config()
function optimize_bridged!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [[128.0, 1.0, -0.0, -32.0, 0.0, 160.0, -0.0, -256.0, 0.0]; zeros(90)],
    )
end
function optimize_cached!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [zeros(90); [128.0, 1.0, -0.0, -32.0, 0.0, 160.0, -0.0, -256.0, 0.0]],
    )
end
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.chebyshev_test(mock, config)
end
