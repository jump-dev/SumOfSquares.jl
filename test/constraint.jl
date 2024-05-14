import SemialgebraicSets
import PolyJuMP

@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    model = SOSModel()
    path = joinpath(
        dirname(dirname(pathof(SumOfSquares))),
        "test",
        "constraint.jl",
    )
    err = ErrorException(
        "At $path:15: `@constraint(model, [1 x; -x 0] in SOSMatrixCone())`: The polynomial matrix constrained to be SOS must be symmetric.",
    )
    @test_throws err @constraint(model, [1 x; -x 0] in SOSMatrixCone())
end

@testset "Printing" begin
    @polyvar x
    model = SOSModel()
    @variable(model, a)
    cref = @constraint(model, a * x^2 >= 1)
    @test sprint(show, MIME"text/plain"(), cref) == "(-1) + (a)x² is SOS"
    @test sprint(show, MIME"text/latex"(), cref) ==
          "\$\$ (-1) + (a)x^{2} \\text{ is SOS} \$\$"
    sdref = @constraint(model, a * x^2 in SDSOSCone())
    @test sprint(show, MIME"text/plain"(), sdref) == "(a)x² is SDSOS"
    @test sprint(show, MIME"text/latex"(), sdref) ==
          "\$\$ (a)x^{2} \\text{ is SDSOS} \$\$"
    dref = @constraint(model, a * x^2 in DSOSCone())
    @test sprint(show, MIME"text/plain"(), dref) == "(a)x² is DSOS"
    @test sprint(show, MIME"text/latex"(), dref) ==
          "\$\$ (a)x^{2} \\text{ is DSOS} \$\$"
    model = Model()
    @variable(model, a)
    for sparsity in [Sparsity.NoPattern(), Sparsity.Variable()]
        cref_fix = @constraint(
            model,
            a * x^2 >= 1,
            SOSCone(),
            domain = (@set x == 1),
            sparsity = sparsity
        )
        @test sprint(show, MIME"text/plain"(), cref_fix) ==
              "(-1) + (a)x² is SOS"
        @test sprint(show, MIME"text/latex"(), cref_fix) ==
              "\$\$ (-1) + (a)x^{2} \\text{ is SOS} \$\$"
        cref_fix = @constraint(
            model,
            a * x^2 >= 1,
            SOSCone(),
            domain = (@set x >= 1),
            sparsity = sparsity
        )
        @test sprint(show, MIME"text/plain"(), cref_fix) ==
              "(-1) + (a)x² is SOS"
        @test sprint(show, MIME"text/latex"(), cref_fix) ==
              "\$\$ (-1) + (a)x^{2} \\text{ is SOS} \$\$"
    end
end

@testset "Bridges with complex numbers" begin
    @polyvar x y
    p = (x + im * y) * (x - im * y)
    model = Model()
    cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
    @constraint(model, p in cone)
    @test SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{ComplexF64} in
          model.bridge_types
    @test SumOfSquares.Bridges.Variable.KernelBridge{ComplexF64} in
          model.bridge_types
end
