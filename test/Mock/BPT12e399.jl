config = MOI.Test.TestConfig()
optimize!_max(mock) = MOIU.mock_optimize!(mock, [ 6.0, 1.0, 9.0, -3.0*√2],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, 1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, √2]]
)
optimize!_min(mock) = MOIU.mock_optimize!(mock, [-6.0, 1.0, 9.0, 3.0*√2],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, -1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, -√2]]
)
mock = bridged_mock(optimize!_max, optimize!_min)
Tests.BPT12e399_test(mock, config)
