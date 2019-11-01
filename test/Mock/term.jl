config = MOI.Test.TestConfig()

@testset "Model" begin
    optimize!(mock) = MOIU.mock_optimize!(mock, [0.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]])
    for mock in mocks(optimize!)
        Tests.sos_term_test(mock, config)
        Tests.sdsos_term_test(mock, config)
        Tests.dsos_term_test(mock, config)
    end
end

# The VectorOfVariables-in-SOSPolynomialSet is bridged by VectorFunctionizeBridge
# since the free variable is bridged. This tests that the GramMatrixAttribute, ...
# are passed by the VectorFunctionizeBridge.
@testset "NoFreeVariable" begin
    optimize!(mock) = MOIU.mock_optimize!(mock, [0.0, 0.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]])
    nofree_mock = bridged_mock(optimize!, model = NoFreeVariable{Float64}())
    for mock in mocks(optimize!)
        Tests.sos_term_test(nofree_mock, config)
        Tests.sdsos_term_test(nofree_mock, config)
        Tests.dsos_term_test(nofree_mock, config)
    end
end
