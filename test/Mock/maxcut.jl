config = MOI.Test.Config()
function optimize1!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        MOI.INFEASIBLE,
        MOI.NO_SOLUTION,
        MOI.INFEASIBILITY_CERTIFICATE,
    )
end
# The test does not check the solution so we just set zeros.
function optimize2!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        zeros(MOI.get(mock, MOI.NumberOfVariables())),
    )
end
for mock in mocks(optimize1!, optimize2!)
    Tests.maxcut_test(mock, config)
end
