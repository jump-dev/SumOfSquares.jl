config = MOI.Test.Config()
# The test does not check the solution except the first variable so we just set zeros.
function optimize_bridged!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-0.184667; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)],
    )
end
function optimize_cached!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1); -0.184667],
    )
end
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sos_quartic_comparison_test(mock, config)
end
function optimize_bridged!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-3.172412; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)],
    )
end
function optimize_cached!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1); -3.172412],
    )
end
for mock in [bridged_mock(optimize_bridged!), cached_mock(optimize_cached!)]
    Tests.sdsos_quartic_comparison_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-11 / 3; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)],
    )
end
for mock in mocks(optimize!)
    Tests.dsos_quartic_comparison_test(mock, config)
end
