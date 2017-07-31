@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    m = SOSModel()
    @test_throws ArgumentError addpolyconstraint!(m, [1 x; -x 0], PSDCone(), BasicSemialgebraicSet())
end
