config = MOI.Test.TestConfig()
optimize!_max(mock) = MOIU.mock_optimize!(mock, [ 6.0, 1.0, 9.0, -3.0*√2],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, 1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, √2]]
)
optimize!_min(mock) = MOIU.mock_optimize!(mock, [-6.0, 1.0, 9.0, 3.0*√2],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, -1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, -√2]]
)
for mock in mocks(optimize!_max, optimize!_min)
    Tests.BPT12e399_rem_test(mock, config)
end
optimize!_max(mock) = MOIU.mock_optimize!(mock, [ 10.0, 4.0, 0.0, 5.0, 0.0, -5.0, 5.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.0, 1, 0, 1, 1]]
)
optimize!_min(mock) = MOIU.mock_optimize!(mock, [-10.0, 4.0, 0.0, 5.0, 0.0,  5.0, 5.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.0, 1, 0, -1, 1]]
)
for mock in mocks(optimize!_max, optimize!_min)
    Tests.BPT12e399_maxdegree_test(mock, config)
end
