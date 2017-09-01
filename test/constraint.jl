@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    m = SOSModel()
    # Throws an Argument Error because it should be symmetric to be an SOS matrix
    @test_throws ArgumentError addpolyconstraint!(m, [1 x; -x 0], PSDCone(), BasicSemialgebraicSet{Int, polynomialtype(x, Int)}())
end
