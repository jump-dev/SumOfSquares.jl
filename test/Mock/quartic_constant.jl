config = MOI.Test.Config()
optimize!(mock) = MOI.Utilities.mock_optimize!(mock, [4.0, 1.0])
for mock in mocks(optimize!)
    Tests.dsos_quartic_constant_test(mock, config)
    Tests.sdsos_quartic_constant_test(mock, config)
    Tests.sos_quartic_constant_test(mock, config)
end
