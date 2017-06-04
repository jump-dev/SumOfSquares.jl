@testset "Creating polynomial with empty MonomialVector" begin
    x = MonomialVector{true}()
    m = SOSModel()
    @inferred createnonnegativepoly(m, Poly{true, :Gram}(x), :Cont)
    @test createnonnegativepoly(m, Poly{true, :Gram}(x), :Cont) == 0
end
