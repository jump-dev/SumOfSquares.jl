config = MOI.Test.Config()
optimize1!(mock) = MOI.Utilities.mock_optimize!(mock, MOI.INFEASIBLE, MOI.NO_SOLUTION, MOI.INFEASIBILITY_CERTIFICATE)
# The test does not check the solution so we just set zeros.
optimize2!(mock) = MOI.Utilities.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
for mock in mocks(optimize1!, optimize2!)
    Tests.motzkin_test(mock, config)
end
