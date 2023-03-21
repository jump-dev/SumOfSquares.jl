config = MOI.Test.Config()
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        MOI.INFEASIBLE,
        MOI.NO_SOLUTION,
        MOI.INFEASIBILITY_CERTIFICATE,
    )
end
for mock in mocks(optimize!)
    Tests.sosdemo5_infeasible_test(mock, config)
end
# The test does not check the solution so we just set zeros.
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        zeros(MOI.get(mock, MOI.NumberOfVariables())),
    )
end
for mock in mocks(optimize!)
    Tests.sosdemo5_feasible_test(mock, config)
end
