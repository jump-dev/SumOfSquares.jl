config = MOI.Test.Config()
function optimize!_max(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [6.0, 9.0, 1.0, -3.0 * √2],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1 / 3, 1, 3]],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) => [[3, 1 / 3, √2]],
    )
end
function optimize!_min(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-6.0, 9.0, 1.0, 3.0 * √2],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1 / 3, -1, 3]],
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone) =>
            [[3, 1 / 3, -√2]],
    )
end
for mock in mocks(optimize!_max, optimize!_min)
    Tests.BPT12e399_rem_test(mock, config)
end
function optimize!_max(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [10.0, 5.0, -5.0, 5.0, 0.0, 0.0, 4.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1, 1, 0.0, 1, 0]],
    )
end
function optimize!_min(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [-10.0, 5.0, 5.0, 5.0, 0.0, 0.0, 4.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1, -1, 0.0, 1, 0]],
    )
end
for mock in [
    bridged_mock(optimize!_max, optimize!_min),
    cached_mock(optimize!_max, optimize!_min),
]
    Tests.BPT12e399_maxdegree_test(mock, config)
end
