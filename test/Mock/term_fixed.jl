config = MOI.Test.Config()
function optimize_bridged!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [1.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]],
    )
end
function optimize_cached!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [0.0, 1.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]],
    )
end
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sos_term_fixed_test(mock, config)
    Tests.sdsos_term_fixed_test(mock, config)
    Tests.dsos_term_fixed_test(mock, config)
end
