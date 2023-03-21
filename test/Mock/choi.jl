config = MOI.Test.Config()
optimize!(mock) = MOI.Utilities.mock_optimize!(
    mock, MOI.INFEASIBLE, MOI.NO_SOLUTION, MOI.INFEASIBILITY_CERTIFICATE)
for mock in mocks(optimize!)
    Tests.choi_test(mock, config)
end
