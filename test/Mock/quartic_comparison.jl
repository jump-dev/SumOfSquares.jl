config = MOI.Test.Config()
# The test does not check the solution except the first variable so we just set zeros.
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-0.184667; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)],
    )
end
for mock in [bridged_mock(optimize!)]
    Tests.sos_quartic_comparison_test(mock, config)
end
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-3.172412; zeros(MOI.get(mock, MOI.NumberOfVariables()) - 1)],
    )
end
for mock in mocks(optimize!)
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
