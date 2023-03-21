config = MOI.Test.Config(atol = 1e-5, rtol = 1e-5)
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0],
    )
end
for mock in mocks(optimize!)
    Tests.dsos_univariate_sum_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(mock, [1.0, 1.0, √2, 1.0, 1.0, -√2])
end
for mock in mocks(optimize!)
    Tests.sdsos_univariate_sum_test(mock, config)
    Tests.sos_univariate_sum_test(mock, config)
end
