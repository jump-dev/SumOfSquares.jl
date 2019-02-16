import SemialgebraicSets
import PolyJuMP

@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    model = SOSModel()
    err = ErrorException("In `@constraint(model, [1 x; -x 0] in SOSMatrixCone())`: The polynomial matrix constrained to be SOS must be symmetric.")
    @test_throws err @constraint(model, [1 x; -x 0] in SOSMatrixCone())
end

@testset "Printing" begin
    @polyvar x
    model = SOSModel()
    @variable(model, a)
    cref = @constraint(model, a * x^2 >= 1)
    @test sprint(show, MIME"text/plain"(), cref) == "(a)x² + (-1) ∈ SOSCone()"
    @test sprint(show, MIME"text/latex"(), cref) == "\$ (a)x^{2} + (-1) \\in SOSCone() \$"
end
