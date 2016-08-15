# Adapted from:
# SOSDEMO2 --- Lyapunov Function Search
# Section 3.2 of SOSTOOLS User's Manual

@polyvar x1 x2 x3

# Constructing the vector field dx/dt = f
f = [-x1^3-x1*x3^2,
    -x2-x1^2*x2,
    -x3+3*x1^2*x3-3*x3/(x3^2+1)]

m = JuMP.Model(solver = SCSSolver(verbose=false))

# The Lyapunov function V(x):
@SOSvariable m V >= 0 [x1^2, x2^2, x3^2]
[prog,V] = sospolyvar(prog, [x1^2, x2^2, x3^2])

@SOSconstraint m V >= x1^2+x2^2+x3^2

# dV/dx*(x3^2+1)*f <= 0
@SOSconstraint (diff(V,x1)*f(1)+diff(V,x2)*f(2)+diff(V,x3)*f(3))*(x3^2+1) <= 0

status = solve(m)

@test status == :Optimal
