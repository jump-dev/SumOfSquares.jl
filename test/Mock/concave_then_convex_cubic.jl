config = MOI.Test.TestConfig()
optimize!(mock) = MOIU.mock_optimize!(mock, [1.0; zeros(5); 6.0; zeros(12); 6.0; zeros(10)],
    (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) => zeros(26),
    (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(7), zeros(7)],
)
for mock in mocks(optimize!)
    Tests.dsos_concave_then_convex_cubic_test(mock, config)
end
optimize!(mock) = MOIU.mock_optimize!(mock, [1.0; zeros(5); 3.0; zeros(3); 3.0; zeros(7); 3.0; zeros(3); 3.0; zeros(5)],
    (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(7), zeros(7)],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [zeros(3) for i in 1:8],
)
for mock in mocks(optimize!)
    Tests.sdsos_concave_then_convex_cubic_test(mock, config)
end
optimize!(mock) = MOIU.mock_optimize!(mock, [1.0; zeros(5); 6.0; zeros(8); 6.0; zeros(6)],
    (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) => zeros(4),
    (MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle) => [zeros(6), zeros(6)],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(7), zeros(7)],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [zeros(3), zeros(3)],
)
for mock in mocks(optimize!)
    Tests.sos_concave_then_convex_cubic_test(mock, config)
end
