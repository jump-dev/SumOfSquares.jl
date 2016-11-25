# Example 3.95 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
facts("Example 3.25") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    @polyvar w x y z
    p = (w^4+1)*(x^4+1)*(y^4+1)*(z^4+1)+2w+3x+4y+5z
    m = Model(solver=solver)
    c = @polyconstraint m p >= 0
    status = solve(m)
    @fact status --> :Optimal
    #s = getslack(c)
    #@show length(s.x)
end; end; end
