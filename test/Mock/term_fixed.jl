config = MOI.Test.Config()
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [1.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]],
    )
end
for mock in mocks(optimize!)
    Tests.sos_term_fixed_test(mock, config)
    Tests.sdsos_term_fixed_test(mock, config)
    Tests.dsos_term_fixed_test(mock, config)
end
