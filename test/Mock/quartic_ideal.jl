config = MOI.Test.Config()
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        MOI.INFEASIBLE,
        MOI.NO_SOLUTION,
        MOI.INFEASIBILITY_CERTIFICATE,
    )
end
@testset "quartic_ideal_rem $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_rem_test(mock, config)
end
@testset "quartic_ideal_2_rem $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_2_rem_test(mock, config)
end
# The test does not check the solution so we just set zeros.
function optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        zeros(MOI.get(mock, MOI.NumberOfVariables())),
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [zeros(3)],
    )
end
@testset "quartic_ideal $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_test(mock, config)
end
@testset "quartic_ideal_4 $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_4_test(mock, config)
end
@testset "quartic_ideal_4_rem $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_4_rem_test(mock, config)
end
