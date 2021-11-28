config = MOI.Test.Config()
optimize!_max_bridged(mock) = MOIU.mock_optimize!(mock, [ 6.0, 1.0, 9.0, -3.0*√2],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, 1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, √2]]
)
optimize!_min_bridged(mock) = MOIU.mock_optimize!(mock, [-6.0, 1.0, 9.0, 3.0*√2],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, -1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, -√2]]
)
optimize!_max_cached(mock) = MOIU.mock_optimize!(mock, [1.0, 9.0, -3.0*√2, 6.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, 1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, √2]]
)
optimize!_min_cached(mock) = MOIU.mock_optimize!(mock, [1.0, 9.0, 3.0*√2, -6.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[3, -1, 1/3]],
    (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1/3, -√2]]
)
for mock in [bridged_mock(optimize!_max_bridged, optimize!_min_bridged), cached_mock(optimize!_max_cached, optimize!_min_cached)]
    Tests.BPT12e399_rem_test(mock, config)
end
optimize!_max_bridged(mock) = MOIU.mock_optimize!(mock, [ 10.0, 4.0, 0.0, 5.0, 0.0, -5.0, 5.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.0, 1, 0, 1, 1]]
)
optimize!_min_bridged(mock) = MOIU.mock_optimize!(mock, [-10.0, 4.0, 0.0, 5.0, 0.0,  5.0, 5.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.0, 1, 0, -1, 1]]
)
optimize!_max_cached(mock) = MOIU.mock_optimize!(mock, [4.0, 0.0, 5.0, 0.0, -5.0, 5.0, 10.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.0, 1, 0, 1, 1]]
)
optimize!_min_cached(mock) = MOIU.mock_optimize!(mock, [4.0, 0.0, 5.0, 0.0,  5.0, 5.0, -10.0],
    (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[0.0, 1, 0, -1, 1]]
)
for mock in [bridged_mock(optimize!_max_bridged, optimize!_min_bridged), cached_mock(optimize!_max_cached, optimize!_min_cached)]
    Tests.BPT12e399_maxdegree_test(mock, config)
end
