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
    Tests.choi_test(mock, config)
end
