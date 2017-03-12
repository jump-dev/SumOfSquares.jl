@testset "Non-symmetric matrix SOS constraint with $solver" for solver in sdp_solvers
    @polyvar x
    m = SOSModel()
    @test_throws ArgumentError addpolynonnegativeconstraint(m, [1 x; -x 0], BasicSemialgebraicSet())
end
