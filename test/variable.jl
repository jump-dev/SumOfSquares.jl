import PolyJuMP

@testset "Creating polynomial with empty MonomialVector" begin
    @polyvar x
    X = emptymonovec(typeof(x^2))
    model = SOSModel()
    v = PolyJuMP.Variable(SOSPoly(X), false, false)
    #@inferred JuMP.add_variable(model, v) # FIXME broken
    @test JuMP.add_variable(model, v) == 0
end
@testset "GramMatrix JuMP.value" begin
    model = JuMP.Model()
    @variable model α
    @variable model β
    @polyvar x y
    q = GramMatrix([α β; β α], [x, y])
    JuMP.fix(α, 2)
    JuMP.fix(β, 3)
    @test_broken JuMP.value(q) == 2x^2 + 2y^2 + 6x*y
end
@testset "Container of GramMatrix" begin
    model = JuMP.Model()
    using DynamicPolynomials
    @polyvar x y
    X = monomials([x, y], 0:2)
    for cone in [SOSPoly(X), SDSOSPoly(X), DSOSPoly(X)]
        p = @variable(model, [1:2], cone)
        @test p[1].x == X
        @test p[2].x == X
        if cone isa SDSOSPoly
            @test eltype(p) <: GramMatrix{JuMP.AffExpr}
        else
            @test eltype(p) <: GramMatrix{JuMP.VariableRef}
        end
    end
end
