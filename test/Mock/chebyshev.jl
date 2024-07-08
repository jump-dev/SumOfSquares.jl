config = MOI.Test.Config()
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [[128.0, 1.0, -0.0, -32.0, 0.0, 160.0, -0.0, -256.0, 0.0]; zeros(90)],
    )
end
for mock in mocks(optimize!)
    Tests.chebyshev_test(mock, config)
end
