# Adapted from:
# SOSDEMO4 --- Matrix Copositivity
# Section 3.4 of SOSTOOLS User's Manual
#
# See also (3.41) of [BPT12].
#
# [BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
# *Semidefinite Optimization and Convex Algebraic Geometry*.
# Society for Industrial and Applied Mathematics, **2012**.


@testset "SOSDEMO4 with $(factory.constructor)" for factory in sdp_factories
    iscsdp(factory) && continue

    # TODO test with CopositiveInner(SOSCone())
    @polyvar x[1:5]

    # Horn matrix
    H = [1 -1  1  1 -1;
        -1  1 -1  1  1;
         1 -1  1 -1  1;
         1  1 -1  1 -1;
        -1  1  1 -1  1]

    xs = vec(x).^2
    xsHxs = dot(xs, H*xs)
    r = sum(xs)

    m0 = SOSModel(factory)
    @constraint m0 xsHxs >= 0
    JuMP.optimize!(m0)
    @test JuMP.dual_status(m0) == MOI.INFEASIBILITY_CERTIFICATE

    m1 = SOSModel(factory)
    @constraint m1 r*xsHxs >= 0
    JuMP.optimize!(m1)
    @test JuMP.primal_status(m1) == MOI.FEASIBLE_POINT
    # Program is feasible. The matrix H is copositive
end
