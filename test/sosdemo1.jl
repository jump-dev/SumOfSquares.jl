# Adapted from:
# SOSDEMO1 --- Sum of Squares Test
# Section 3.1 of SOSTOOLS User's Manual

@polyvar x y

m = JuMP.Model(solver = SCSSolver(verbose=false))

p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4

soscon = @SOSconstraint m p >= 0

status = solve(m)

@test status == :Optimal

# FIXME what should those 3 values be ?
q1 = sqrt(2)*x^2 + sqrt(2)/2*x*y + -0.6692227671373556*y^2
q2 = 0.6263451365977047*x*y + 0.7556127920112614*y^2
q3 = 2*y^2
expected = SOSDecomposition([q1,q2,q3])
@test isapprox(expected, SOSDecomposition(getslack(soscon)); rtol=3e-3)

M = JuMP.Model(solver = SCSSolver(verbose=false))

p = 4*x^4*y^6 + x^2 - x*y^2 + y^2

soscon = @SOSconstraint M p >= 0

status = solve(M)

@test status == :Optimal

# p should be MatPolynomial([1, 0, -1/2, 0, -1, 1, 0, -2/3, 0, 4/3, 0, 0, 2, 0, 4], [y, x, x*y, x*y^2, x^2*y^3])
