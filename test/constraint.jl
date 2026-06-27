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
    @test sprint(show, MIME"text/plain"(), cref) == "(-1)·1 + (a)·x² is SOS"
    @test sprint(show, MIME"text/latex"(), cref) ==
          "\$\$ (-1) \\cdot 1 + (a) \\cdot x^{2} \\text{ is SOS} \$\$"
    sdref = @constraint(model, a * x^2 in SDSOSCone())
    @test sprint(show, MIME"text/plain"(), sdref) == "(a)·x² is SDSOS"
    @test sprint(show, MIME"text/latex"(), sdref) ==
          "\$\$ (a) \\cdot x^{2} \\text{ is SDSOS} \$\$"
    dref = @constraint(model, a * x^2 in DSOSCone())
    @test sprint(show, MIME"text/plain"(), dref) == "(a)·x² is DSOS"
    @test sprint(show, MIME"text/latex"(), dref) ==
          "\$\$ (a) \\cdot x^{2} \\text{ is DSOS} \$\$"
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
              "(-1)·1 + (a)·x² is SOS"
        @test sprint(show, MIME"text/latex"(), cref_fix) ==
              "\$\$ (-1) \\cdot 1 + (a) \\cdot x^{2} \\text{ is SOS} \$\$"
        cref_fix = @constraint(
            model,
            a * x^2 >= 1,
            SOSCone(),
            domain = (@set x >= 1),
            sparsity = sparsity
        )
        @test sprint(show, MIME"text/plain"(), cref_fix) ==
              "(-1)·1 + (a)·x² is SOS"
        @test sprint(show, MIME"text/latex"(), cref_fix) ==
              "\$\$ (-1) \\cdot 1 + (a) \\cdot x^{2} \\text{ is SOS} \$\$"
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

@testset "Generic value type `$T`" for T in [BigFloat, Float32]
    @polyvar x y
    model = GenericModel{T}()
    # The coefficients of the polynomial are converted to the value type `T` of
    # the model, like JuMP does for affine and quadratic constraints. Here, the
    # constant `one(T)` is a `Real`, not a JuMP scalar, and used to leak `Float64`.
    @constraint(model, one(T) + x^2 in SOSCone())
    F = JuMP.GenericAffExpr{T,JuMP.GenericVariableRef{T}}
    @test any(list_of_constraint_types(model)) do (Fi, Si)
        return Fi == Vector{F} && Si <: SumOfSquares.SOSPolynomialSet
    end
    @test SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{T} in
          model.bridge_types
    @test SumOfSquares.Bridges.Variable.KernelBridge{T} in model.bridge_types
    # The `domain` polynomials used to leak `Float64` into the multiplier bridge.
    @constraint(model, one(T) - x^2 in SOSCone(), domain = @set y >= one(T))
    @test SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge{T} in
          model.bridge_types
end
