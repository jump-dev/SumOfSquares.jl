# Example 3.77 and 3.79 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
facts("Example 3.77 and 3.79") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    @polyvar x
    P = [x^2-2x+2 x; x x^2]
    # Example 3.77
    m = Model(solver=solver)
    @polyconstraint m P >= 0
    status = solve(m)
    @fact status --> :Optimal
    # Example 3.79
    @polyvar y[1:2]
    M = Model(solver=solver)
    @polyconstraint M dot(y, P*y) >= 0
    status = solve(M)
    @fact status --> :Optimal
end; end; end
