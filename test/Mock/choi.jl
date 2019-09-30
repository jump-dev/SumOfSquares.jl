config = MOI.Test.TestConfig()
optimize!(mock) = MOIU.mock_optimize!(mock, MOI.INFEASIBLE, tuple(),
                                      MOI.INFEASIBILITY_CERTIFICATE)
for mock in mocks(optimize!)
    Tests.choi_test(mock, config)
end
