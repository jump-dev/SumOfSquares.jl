config = MOI.Test.Config()
# The test does not check the solution so we just set zeros.
optimize!(mock) = MOI.Utilities.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
for mock in mocks(optimize!)
    Tests.simple_matrix_test(mock, config)
end
