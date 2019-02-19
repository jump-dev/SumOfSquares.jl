# Test that the macro call `m` throws an error exception during pre-compilation
macro test_macro_throws(errortype, m)
    # See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533/5878
    :(@test_throws $errortype try @eval $m catch err; throw(err.error) end)
end

using JuMP

function bridged_mock(mock_optimize!::Function;
                      model = JuMP._MOIModel{Float64}())
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!)
    return bridged
end

function test_noc(bridged_mock, F, S, n)
    @test MOI.get(bridged_mock, MOI.NumberOfConstraints{F, S}()) == n
    @test length(MOI.get(bridged_mock, MOI.ListOfConstraintIndices{F, S}())) == n
    @test ((F, S) in MOI.get(bridged_mock, MOI.ListOfConstraints())) == !iszero(n)
end

# Test deletion of bridge
function test_delete_bridge(m::MOI.Bridges.AbstractBridgeOptimizer,
                            ci::MOI.ConstraintIndex{F, S}, nvars::Int,
                            nocs::Tuple; last_bridge = true) where {F, S}
    @test MOI.get(m, MOI.NumberOfVariables()) == nvars
    test_noc(m, F, S, 1)
    for noc in nocs
        test_noc(m, noc...)
    end
    @test MOI.is_valid(m, ci)
    MOI.delete(m, ci)
    @test_throws MOI.InvalidIndex{typeof(ci)} MOI.delete(m, ci)
    try
        MOI.delete(m, ci)
    catch err
        @test err.index == ci
    end
    @test !MOI.is_valid(m, ci)
    if last_bridge
        @test isempty(m.bridges)
    end
    test_noc(m, F, S, 0)
    # As the bridge has been removed, if the constraints it has created where not removed, it wouldn't be there to decrease this counter anymore
    @test MOI.get(m, MOI.NumberOfVariables()) == nvars
    for noc in nocs
        test_noc(m, noc...)
    end
end
