@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    MOI.empty!(solver)
    m = SOSModel()
    # Throws an Argument Error because it should be symmetric to be an SOS matrix
    @test_throws ArgumentError PolyJuMP.addpolyconstraint!(m, [1 x; -x 0], SOSMatrixCone(), BasicSemialgebraicSet{Int, polynomialtype(x, Int)}(), PolyJuMP.MonomialBasis)
end
