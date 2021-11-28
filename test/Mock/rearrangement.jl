config = MOI.Test.Config(atol=1e-5, rtol=1e-5)
x = [1.0, -1.0, 1.0, 0.0, 0.0, 0.0]
function optimize!(mock)
    MOIU.mock_optimize!(mock, [zeros(MOI.get(mock, MOI.NumberOfVariables()) - 12); x; x],
        (MOI.VectorAffineFunction{Float64},MOI.Zeros) => [zeros(16)])
end
for mock in mocks(optimize!)
    Tests.rearrangement_test(mock, config)
end
