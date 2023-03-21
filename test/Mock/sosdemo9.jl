config = MOI.Test.Config()
# The test does not check the solution so we just set zeros.
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        zeros(MOI.get(mock, MOI.NumberOfVariables())),
    )
end
for mock in mocks(optimize!)
    Tests.sosdemo9_test(mock, config)
end
