config = MOI.Test.TestConfig()
optimize!(mock) = MOIU.mock_optimize!(mock, MOI.INFEASIBLE, tuple(),
                                      MOI.INFEASIBILITY_CERTIFICATE)
mock = bridged_mock(optimize!)
Tests.choi_test(mock, config)
