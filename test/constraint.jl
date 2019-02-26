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
    @test sprint(show, MIME"text/plain"(), cref) == "(a)x² + (-1) is SOS"
    @test sprint(show, MIME"text/latex"(), cref) == "\$ (a)x^{2} + (-1) \\text{ is SOS} \$"
    sdref = @constraint(model, a * x^2 in SDSOSCone())
    @test sprint(show, MIME"text/plain"(), sdref) == "(a)x² is SDSOS"
    @test sprint(show, MIME"text/latex"(), sdref) == "\$ (a)x^{2} \\text{ is SDSOS} \$"
    dref = @constraint(model, a * x^2 in DSOSCone())
    @test sprint(show, MIME"text/plain"(), dref) == "(a)x² is DSOS"
    @test sprint(show, MIME"text/latex"(), dref) == "\$ (a)x^{2} \\text{ is DSOS} \$"
end
