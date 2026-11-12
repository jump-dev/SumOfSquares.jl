using Test, SumOfSquares, DynamicPolynomials
import MathOptInterface as MOI
import SemialgebraicSets
import PolyJuMP
import MultivariateBases as MB
import MultivariatePolynomials as MP
import StarAlgebras as SA

@testset "Slack bridge registration" begin
    for S in (
        SumOfSquares.WeightedSOSCone,
        SumOfSquares.ScaledDiagonallyDominantConeTriangle,
    )
        for F in (MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64})
            for (BT, T) in PolyJuMP.bridges(F, S)
                @test MOI.supports_constraint(BT{T}, F, S)
            end
        end
    end
end

@testset "Certificate basis index mapping" begin
    @polyvar x
    attr = SumOfSquares.CertificateBasis()
    index_map = MOI.Utilities.IndexMap()
    for B in (Monomial, ScaledMonomial, Chebyshev)
        basis = SubBasis{B}([x^0, x])
        for value in (basis, [basis], [[basis]])
            @test MOI.Utilities.map_indices(index_map, attr, value) === value
            @test MOI.Utilities.map_indices(identity, attr, value) === value
        end
    end
end

@testset "Non-symmetric matrix SOS constraint" begin
    @polyvar x
    model = SOSModel()
    try
        @constraint(model, [1 x; -x 0] in SOSMatrixCone())
        @assert false
    catch err
        @test occursin(
            "The polynomial matrix constrained to be SOS must be symmetric.",
            sprint(showerror, err),
        )
    end
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
    @test SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge{
        T,
    } in model.bridge_types
end

@testset "Domain promotion preserves $B" for B in (
    MP.Monomial,
    MB.ScaledMonomial,
    MB.Chebyshev,
)
    @polyvar x y
    basis = MB.FullBasis{B}([y])
    p = MB.convert_basis(basis, 1.0 + y^2)
    gram = GramMatrix([1.0 0.0; 0.0 1.0], MB.SubBasis{B}([y^0, y]))
    decomposition = SOSDecomposition([one(p), MB.convert_basis(basis, 1.0y)])
    for operand in (p, gram, decomposition)
        expected = MP.polynomial(operand)
        for domain in (FullSpace(), @set(x^2 == 1), @set(x >= 0))
            promoted_domain, promoted = SA.promote_bases(domain, operand)
            vars = domain isa FullSpace ? [y] : [x, y]
            @test promoted == expected
            @test eltype(SA.basis(promoted)) <: MP.BasisElement{B}
            @test MP.variables(promoted) == vars
            @test MP.variables(expected) == [y]
            if !(domain isa FullSpace)
                @test MP.variables(promoted_domain) == vars
            end
            # The semialgebraic bridge still requires monomial-basis inputs.
            # Promotion itself supports every basis family tested here.
            if domain isa BasicSemialgebraicSet && B !== MP.Monomial
                continue
            end
            model = SOSModel()
            cref = @constraint(model, operand in SOSCone(), domain = domain)
            obj = JuMP.constraint_object(cref)
            output = JuMP.reshape_vector(obj.func, JuMP.shape(obj))
            @test output == expected
            @test typeof(SA.object(parent(output))) ===
                  typeof(SA.object(parent(expected)))
            @test MP.variables(output) == vars
            @test MP.variables(expected) == [y]
            if !(domain isa FullSpace)
                @test MP.variables(obj.set.domain) == vars
            end
        end
    end
end
