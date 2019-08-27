config = MOI.Test.TestConfig()
# The test does not check the solution so we just set zeros.
optimize!(mock) = MOIU.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
mock = bridged_mock(optimize!)
Tests.sosdemo9_test(mock, config)
