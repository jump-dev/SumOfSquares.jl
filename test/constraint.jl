import SemialgebraicSets
import PolyJuMP

@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    model = SOSModel()
    err = ErrorException("In `@constraint(model, [1 x; -x 0] in SOSMatrixCone())`: The polynomial matrix constrained to be SOS must be symmetric.")
    @test_throws err @constraint(model, [1 x; -x 0] in SOSMatrixCone())
end
