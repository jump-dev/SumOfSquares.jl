config = MOI.Test.TestConfig()
optimize1!(mock) = MOIU.mock_optimize!(mock, MOI.INFEASIBLE, tuple(), MOI.INFEASIBILITY_CERTIFICATE)
# The test does not check the solution so we just set zeros.
optimize2!(mock) = MOIU.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
mock = bridged_mock(optimize1!, optimize2!)
Tests.maxcut_test(mock, config)
