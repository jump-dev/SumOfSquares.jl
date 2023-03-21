# TODO Make it work with `optimize!`
config = MOI.Test.Config(exclude=Any[MOI.optimize!])
_optimize!(mock) = MOI.Utilities.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
for mock in mocks(_optimize!)
    Tests.dsos_options_pricing_test(mock, config)
    Tests.sdsos_options_pricing_test(mock, config)
    Tests.sos_options_pricing_test(mock, config)
end
