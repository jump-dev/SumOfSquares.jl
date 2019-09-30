config = MOI.Test.TestConfig()
optimize!(mock) = MOIU.mock_optimize!(mock, [[128.0, 0.0, -256.0, 0.0, 160.0, 0.0, -32.0, 0.0, 1.0]; zeros(90)])
for mock in mocks(optimize!)
    Tests.chebyshev_test(mock, config)
end
