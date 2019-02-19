include("Tests/Tests.jl")

using Test, JuMP, GLPK
factory = with_optimizer(GLPK.Optimizer)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
