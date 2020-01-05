optimize!(mock) = MOIU.mock_optimize!(mock, MOI.INFEASIBLE, tuple(), MOI.INFEASIBILITY_CERTIFICATE)
@testset "quartic_ideal_rem $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_rem_test(mock, config)
end
@testset "quartic_ideal_2_rem $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_2_rem_test(mock, config)
end
# The test does not check the solution so we just set zeros.
optimize!(mock) = MOIU.mock_optimize!(mock, zeros(MOI.get(mock, MOI.NumberOfVariables())))
@testset "quartic_ideal $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_test(mock, config)
end
@testset "quartic_ideal_4 $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_4_test(mock, config)
end
@testset "quartic_ideal_4_rem $(typeof(mock))" for mock in mocks(optimize!)
    Tests.quartic_ideal_4_rem_test(mock, config)
end
