using JuMP
const MOIT = MOI.Test

function bridged_mock(mock_optimize!::Function...;
                      model = MOI.Utilities.Model{Float64}())
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end
