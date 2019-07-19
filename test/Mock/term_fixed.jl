config = MOI.Test.TestConfig()
optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 0.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]])
mock = bridged_mock(optimize!)
Tests.sos_term_fixed_test(mock, config)
Tests.sdsos_term_fixed_test(mock, config)
Tests.dsos_term_fixed_test(mock, config)
