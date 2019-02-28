using Test

# Test that the macro call `m` throws an error exception during pre-compilation
macro test_macro_throws(errortype, m)
    # See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533/5878
    :(@test_throws $errortype try @eval $m catch err; throw(err.error) end)
end

using JuMP
const MOIT = MOI.Test

function bridged_mock(mock_optimize!::Function...;
                      model = JuMP._MOIModel{Float64}())
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end
