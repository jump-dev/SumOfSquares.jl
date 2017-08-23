@testset "Creating polynomial with empty MonomialVector" begin
    @polyvar x
    X = emptymonovec(typeof(x^2))
    m = SOSModel()
    @inferred createpoly(m, Poly{true, :Gram}(X), :Cont)
    @test createpoly(m, Poly{true, :Gram}(X), :Cont) == 0
end
@testset "MatPolynomial getvalue" begin
    m = JuMP.Model()
    @variable m α
    @variable m β
    @polyvar x y
    q = MatPolynomial([α β; β α], [x, y])
    JuMP.fix(α, 2)
    JuMP.fix(β, 3)
    @test getvalue(q) == 2x^2 + 2y^2 + 6x*y
end
