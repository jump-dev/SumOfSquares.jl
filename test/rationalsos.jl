# Example 3.45 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
facts("Example 3.45") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    @polyvar x y z w
    p = 2x^4 + x^2*y^2 + y^4 - 4x^2*z - 4x*y*z - 2y^2*w + y^2 - 2y*z + 8z^2 - 2z*w + 2w^2
    m = Model(solver=solver)
    @polyconstraint m p >= 0
    status = solve(m)
    @fact status --> :Optimal
end; end; end
