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
    m = JuMP.Model()
    @variable m α
    @variable m β
    @polyvar x y
    q = GramMatrix([α β; β α], [x, y])
    JuMP.fix(α, 2)
    JuMP.fix(β, 3)
    @test_broken JuMP.value(q) == 2x^2 + 2y^2 + 6x*y
end
