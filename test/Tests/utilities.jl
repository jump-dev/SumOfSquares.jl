using Test, JuMP
const MOIT = MOI.Test
const MOIB = MOI.Bridges

function _model(optimizer::MOI.AbstractOptimizer)
    MOI.empty!(optimizer)
    return direct_model(optimizer)
end

function _model(factory::OptimizerFactory)
    return Model(factory)
end

#"""
#    @test_suite setname subsets
#
#Defines a function `setname_test(model, config, exclude)` that runs the tests
#defined in the dictionary `setname_tests` with the model `model` and config
#`config` except the tests whose dictionary key is in `exclude`. If `subsets` is
#`true` then each test runs in fact multiple tests hence the `exclude` argument
#is passed as it can also contains test to be excluded from these subsets of
#tests.
#"""
macro test_suite(setname, subsets=false)
    testname = Symbol(string(setname) * "_test")
    testdict = Symbol(string(testname) * "s")
    if subsets
        runtest = :( f(model, config, exclude) )
    else
        runtest = :( f(model, config) )
    end
    esc(:(
      function $testname(model::Union{$MOI.ModelLike, OptimizerFactory},
                         config::$MOI.Test.TestConfig,
                         exclude::Vector{String} = String[])
            for (name,f) in $testdict
                if name in exclude
                    continue
                end
                @testset "$name" begin
                    $runtest
                end
            end
        end
    ))
end

function test_noc(model, F, S, n)
    @test MOI.get(model, MOI.NumberOfConstraints{F, S}()) == n
    @test length(MOI.get(model, MOI.ListOfConstraintIndices{F, S}())) == n
    @test ((F, S) in MOI.get(model, MOI.ListOfConstraints())) == !iszero(n)
end

# Test deletion of bridge
function test_delete_bridge(model::Model,
                            cref::ConstraintRef{Model,
                                                MOI.ConstraintIndex{F, S}},
                            nvars::Int, nocs::Tuple;
                            last_bridge = true) where {F, S}
    @test num_variables(model) == nvars
    test_noc(model, F, S, 1)
    for noc in nocs
        test_noc(model, noc...)
    end
    @test is_valid(model, cref)
    delete(model, cref)
    @test_throws MOI.InvalidIndex(index(cref)) delete(model, cref)
    @test !is_valid(model, cref)
    test_noc(model, F, S, 0)
    # As the bridge has been removed, if the constraints it has created where not removed, it wouldn't be there to decrease this counter anymore
    @test num_variables(model) == nvars
    for noc in nocs
        test_noc(model, noc...)
    end
end

# Utilities for building the mock `optimize!` from the solution of a solver
_inner(model::MOIU.CachingOptimizer) = _inner(model.optimizer)
_inner(model::MOIB.LazyBridgeOptimizer) = model.model
# Variables primal values for inner bridged model
function print_value(v, atol)
    i = round(v)
    if isapprox(v, i, atol=atol)
        print(float(i))
    else
        print(v)
    end
end
function inner_variable_value(model, atol=1e-4)
    inner = _inner(backend(model))
    values = MOI.get(inner, MOI.VariablePrimal(),
                     MOI.get(inner, MOI.ListOfVariableIndices()))
    print("optimize!(mock) = MOIU.mock_optimize!(mock, [")
    for (i, v) in enumerate(values)
        if i > 1
            print(", ")
        end
        print_value(v, atol)
    end
    println("])")
    return values
end
# Constraint dual values for inner bridged model
