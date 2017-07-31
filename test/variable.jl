@testset "Creating polynomial with empty MonomialVector" begin
    x = MonomialVector{true}()
    m = SOSModel()
    @inferred createpoly(m, Poly{true, :Gram}(x), :Cont)
    @test createpoly(m, Poly{true, :Gram}(x), :Cont) == 0
end
