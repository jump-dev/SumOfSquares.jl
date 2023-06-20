using JuMP

MOI.Utilities.@model(
    NoFreeVariable,
    (),
    (MOI.EqualTo, MOI.LessThan, MOI.GreaterThan),
    (
        MOI.Nonnegatives,
        MOI.Nonpositives,
        MOI.Zeros,
        MOI.RotatedSecondOrderCone,
        MOI.PositiveSemidefiniteConeTriangle,
    ),
    (),
    (),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,)
)
# No free variables to make sure variable bridges are used to increase coverage
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Reals},
)
    return false
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Nonnegatives},
)
    return true
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Nonpositives},
)
    return true
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Zeros},
)
    return true
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.RotatedSecondOrderCone},
)
    return true
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
)
    return true
end

function bridged_mock(
    mock_optimize!::Function...;
    model = MOI.Utilities.Model{Float64}(),
)
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end

function cached_mock(
    args...;
    cache = () ->
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    kws...,
)
    # We want the MOI backend of the JuMP model to be in `MOI.Utilities.EMPTY_OPTIMIZER`
    # mode so that it's copied at `optimize!` to test that `copy_to` works.
    # If we just return `cached`, it will be emptied in `_model` and the state
    # will be `MOI.Utilities.ATTACHED_OPTIMIZER` which is not what we want. For this
    # reason we return a `JuMP.OptimizerFactory` which returns `cached` instead.
    return (
        () -> begin
            cached = MOI.Utilities.CachingOptimizer(
                cache(),
                MOI.Utilities.AUTOMATIC,
            )
            optimizer = bridged_mock(args...; kws...)
            MOI.Utilities.reset_optimizer(cached, optimizer)
            return cached
        end
    )
end

function mocks(args...; kws...)
    return [bridged_mock(args...; kws...), cached_mock(args...; kws...)]
end

function cheat(test, optimizer_constructor, args...)
    # We add a first layer of bridges that will only apply the SumOfSquares bridges
    # since the `JuMP.add_bridge`s will only apply to the first layer.
    # That way, we know that after this first layer, it will be the same as the mock model.
    # We add a `CachingOptimizer` inbetween as two layers of bridges might mess up with bridged variable
    # indices as they both interpret negative indices as bridged variables that they have bridged
    function constructor()
        model = MOI.Utilities.CachingOptimizer(
            MOI.Utilities.Model{Float64}(),
            MOI.instantiate(optimizer_constructor, with_bridge_type = Float64),
        )
        @show MOI.Utilities.state(model)
        # We attach to make sure the variables are added in the same order
        # But the bridge costs attributes are still passed along so the order might still differ
        MOI.Utilities.attach_optimizer(model)
        @show MOI.Utilities.state(model)
        return model
    end
    model = test(constructor, args...)
    return Tests.inner_variable_value(model)
end
