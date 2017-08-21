@testset "Creating polynomial with empty MonomialVector" begin
    x = MonomialVector{true}()
    m = SOSModel()
    @inferred createpoly(m, Poly{true, :Gram}(x), :Cont)
    @test createpoly(m, Poly{true, :Gram}(x), :Cont) == 0
end
@testset "MatPolynomial getvalue" begin
    q = MatPolynomial([α β; β α], [x, y])
    @test getvalue(q) == 2x^2 + 2y^2 + 6x*y
end
