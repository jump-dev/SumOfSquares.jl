@testset "Creating polynomial with empty MonomialVector" begin
    x = MonomialVector{true}()
    m = SOSModel()
    @inferred createnonnegativepoly(m, :Gram, x)
    @test createnonnegativepoly(m, :Gram, x) == 0
end
