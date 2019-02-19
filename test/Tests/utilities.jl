using JuMP

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
