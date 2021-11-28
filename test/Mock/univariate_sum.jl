config = MOI.Test.Config(atol=1e-5, rtol=1e-5)
optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0])
for mock in mocks(optimize!)
    Tests.dsos_univariate_sum_test(mock, config)
end
optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0, √2, 1.0, 1.0, -√2])
for mock in mocks(optimize!)
    Tests.sdsos_univariate_sum_test(mock, config)
    Tests.sos_univariate_sum_test(mock, config)
end
