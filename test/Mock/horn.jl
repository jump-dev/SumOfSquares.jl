config = MOI.Test.Config()
optimize!_inf(mock) = MOI.Utilities.mock_optimize!(mock, MOI.INFEASIBLE)
optimize!(mock) = MOI.Utilities.mock_optimize!(mock,
    [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0,
     0.0, 0.0, 0.0, 0.0, -0.0, 0.0,
     1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0,
     4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, -0.0, 0.0, -0.0, 0.0, 0.0, -0.0, 0.0])
@testset "sos $(typeof(mock))" for mock in mocks(optimize!_inf, optimize!_inf, optimize!)
    Tests.sos_horn_test(mock, config)
end
@testset "dsos $(typeof(mock))" for mock in mocks(optimize!_inf)
    Tests.dsos_horn_test(mock, config)
end
@testset "sdsos $(typeof(mock))" for mock in mocks(optimize!_inf)
    Tests.sdsos_horn_test(mock, config)
end
