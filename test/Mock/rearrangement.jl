config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5)
x = [1.0, -1.0, 1.0, 0.0, 0.0, 0.0]
optimize!(mock) = MOIU.mock_optimize!(mock, [zeros(MOI.get(mock, MOI.NumberOfVariables()) - 12); x; x])
for mock in mocks(optimize!)
    Tests.rearrangement_test(mock, config)
end
