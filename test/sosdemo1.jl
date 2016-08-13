# Adapted from:
# SOSDEMO1 --- Sum of Squares Test
# Section 3.1 of SOSTOOLS User's Manual

@polyvar x1 x2

m = JuMP.Model(solver = SCSSolver(verbose=false))

p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4

@SOSconstraint m p >= 0

status = solve(m)

@test status == :Optimal
