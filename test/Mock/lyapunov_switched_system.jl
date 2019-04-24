config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5)
α = 0.7291971148804172
γ = √2 + 2e-1
β = (γ^2 - 2) / γ^2 * (1 + α)
# We have the solution
# P_0 = [α 0
#        0 α]
# P   = [1+α 0
#        0   1+α]
# and the slack of the LMI's
# P - A1'P*A*1 = [β   0
#                 0   1+α]
# P - A2'P*A*2 = [1+α 0
#                 0   β]
optimize!(mock) = MOIU.mock_optimize!(mock, [α, 0.0, α, 1 + α, 0.0, β, 1 + α, β, 0.0])
mock = bridged_mock(optimize!)
Tests.quadratic_feasible_lyapunov_switched_system_test(mock, config)
# TODO quadratic infeasible and quartic
