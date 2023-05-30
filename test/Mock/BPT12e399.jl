config = MOI.Test.Config()
function optimize!_max_bridged(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [6.0, 9.0, 1.0, -3.0 * √2],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1 / 3, 1, 3]],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1 / 3, √2]],
    )
end
function optimize!_min_bridged(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-6.0, 9.0, 1.0, 3.0 * √2],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1 / 3, -1, 3]],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [[3, 1 / 3, -√2]],
    )
end
function optimize!_max_cached(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [9.0, 1.0, -3.0 * √2, 6.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1 / 3, 1, 3]],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[1 / 3, 3, √2]],
    )
end
function optimize!_min_cached(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [9.0, 1.0, 3.0 * √2, -6.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1 / 3, -1, 3]],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [[1 / 3, 3, -√2]],
    )
end
for mock in [
    bridged_mock(optimize!_max_bridged, optimize!_min_bridged),
    cached_mock(optimize!_max_cached, optimize!_min_cached),
]
    Tests.BPT12e399_rem_test(mock, config)
end
function optimize!_max_bridged(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [10.0, 5.0, -5.0, 5.0, 0.0, 0.0, 4.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1, 1, 0.0, 1, 0]],
    )
end
function optimize!_min_bridged(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-10.0, 5.0, 5.0, 5.0, 0.0, 0.0, 4.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1, -1, 0.0, 1, 0]],
    )
end
function optimize!_max_cached(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [5.0, -5.0, 5.0, 0.0, 0.0, 4.0, 10.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1, 1, 0.0, 1, 0]],
    )
end
function optimize!_min_cached(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [5.0, 5.0, 5.0, 0.0, 0.0, 4.0, -10.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1, -1, 0.0, 1, 0]],
    )
end
for mock in [
    bridged_mock(optimize!_max_bridged, optimize!_min_bridged),
    cached_mock(optimize!_max_cached, optimize!_min_cached),
]
    Tests.BPT12e399_maxdegree_test(mock, config)
end
