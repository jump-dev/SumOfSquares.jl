using JuMP
const MOIT = MOI.Test

MOIU.@model(NoFreeVariable,
            (), (MOI.EqualTo, MOI.LessThan, MOI.GreaterThan), (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros, MOI.RotatedSecondOrderCone, MOI.PositiveSemidefiniteConeTriangle), (),
            (), (MOI.ScalarAffineFunction,), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,))
# No free variables to make sure variable bridges are used to increase coverage
MOI.supports_constraint(::NoFreeVariable, ::Type{MOI.VectorOfVariables}, ::Type{MOI.Reals}) = false

function bridged_mock(mock_optimize!::Function...;
                      model = MOI.Utilities.Model{Float64}())
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end

function cached_mock(
    args...; cache = () -> MOIU.UniversalFallback(MOIU.Model{Float64}()), kws...)
    # We want the MOI backend of the JuMP model to be in `MOIU.EMPTY_OPTIMIZER`
    # mode so that it's copied at `optimize!` to test that `copy_to` works.
    # If we just return `cached`, it will be emptied in `_model` and the state
    # will be `MOIU.ATTACHED_OPTIMIZER` which is not what we want. For this
    # reason we return a `JuMP.OptimizerFactory` which returns `cached` instead.
    with_optimizer(() -> begin
        cached = MOIU.CachingOptimizer(cache(), MOIU.AUTOMATIC)
        optimizer = bridged_mock(args...; kws...)
        MOIU.reset_optimizer(cached, optimizer)
        return cached
    end)
end

function mocks(args...; kws...)
    return [bridged_mock(args...; kws...), cached_mock(args...; kws...)]
end
