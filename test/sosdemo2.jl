# Adapted from:
# SOSDEMO2 --- Lyapunov Function Search
# Section 3.2 of SOSTOOLS User's Manual
using Calculus

@polyvar x1 x2 x3

# Constructing the vector field dx/dt = f
f = [-x1^3-x1*x3^2,
    -x2-x1^2*x2,
    -x3+3*x1^2*x3-3*x3/(x3^2+1)]

m = JuMP.Model(solver = SCSSolver(verbose=false))

# The Lyapunov function V(x):
@SOSvariable m V >= 0 [x1^2, x2^2, x3^2]

@SOSconstraint m V >= x1^2+x2^2+x3^2

# dV/dx*(x3^2+1)*f <= 0
x = [x1, x2, x3]
P = dot(differentiate(V, x), f)*(x3^2+1)
@SOSconstraint m P <= 0

status = solve(m)

@test status == :Optimal

# SOSTools doc:
# expected = 3.0922 * x1^2 + 2.2885 * x2^2 + x3^2
# What I get with SCS
expected = 3.937089849725135x1^2 + 1.089493395427901x2^2 + 1.1692659364561195x3^2
obtained = VecPolynomial(getvalue(V))
@test isapprox(obtained, expected)
