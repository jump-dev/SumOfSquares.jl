config = MOI.Test.TestConfig()
# The test does not check the solution so we just set zeros.
optimize!(mock) = MOIU.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
for mock in mocks(optimize!)
    Tests.sosdemo10_test(mock, config)
end
