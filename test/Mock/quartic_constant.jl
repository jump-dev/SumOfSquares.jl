config = MOI.Test.Config()
optimize_bridged!(mock) = MOI.Utilities.mock_optimize!(mock, [4.0, 1.0])
optimize_cached!(mock) = MOI.Utilities.mock_optimize!(mock, [1.0, 4.0])
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.dsos_quartic_constant_test(mock, config)
    Tests.sdsos_quartic_constant_test(mock, config)
    Tests.sos_quartic_constant_test(mock, config)
end
