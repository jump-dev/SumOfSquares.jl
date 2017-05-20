@testset "Creating polynomial with empty MonomialVector" begin
    x = MonomialVector{true}()
    m = SOSModel()
    @inferred createnonnegativepoly(m, :Gram, x, :Cont)
    @test createnonnegativepoly(m, :Gram, x, :Cont) == 0
end
