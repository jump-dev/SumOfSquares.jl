config = MOI.Test.Config()
# The test does not check the solution except the first variable so we just set zeros.
optimize_bridged!(mock) = MOIU.mock_optimize!(mock, [-0.184667; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)])
optimize_cached!(mock) = MOIU.mock_optimize!(mock, [zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1); -0.184667])
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sos_quartic_comparison_test(mock, config)
end
optimize_bridged!(mock) = MOIU.mock_optimize!(mock, [-3.172412; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)])
optimize_cached!(mock) = MOIU.mock_optimize!(mock, [zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1); -3.172412])
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sdsos_quartic_comparison_test(mock, config)
end
optimize!(mock) = MOIU.mock_optimize!(mock, [-11/3; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)])
for mock in mocks(optimize!)
    Tests.dsos_quartic_comparison_test(mock, config)
end
