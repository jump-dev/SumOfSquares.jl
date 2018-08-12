@testset "Choi with $(factory.constructor)" for factory in sdp_factories
    @polyvar x y z

    m = SOSModel(factory)

    C = [x^2+2y^2 -x*y -x*z;
         -x*y y^2+2z^2 -y*z;
         -x*z -y*z z^2+2x^2]

    @SDconstraint m C >= 0

    JuMP.optimize(m)
    @test JuMP.dualstatus(m) == MOI.InfeasibilityCertificate
end
