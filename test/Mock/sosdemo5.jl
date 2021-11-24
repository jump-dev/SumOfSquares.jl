config = MOI.Test.Config()
optimize!(mock) = MOIU.mock_optimize!(mock, MOI.INFEASIBLE, tuple(), MOI.INFEASIBILITY_CERTIFICATE)
for mock in mocks(optimize!)
    Tests.sosdemo5_infeasible_test(mock, config)
end
# The test does not check the solution so we just set zeros.
optimize!(mock) = MOIU.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
for mock in mocks(optimize!)
    Tests.sosdemo5_feasible_test(mock, config)
end
