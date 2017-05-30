@testset "Monomial selection for certificate" begin
    @polyvar x y
    @test_throws ErrorException getmonomialsforcertificate([x*y, y^2], :Sparse)
    @test getmonomialsforcertificate([x*y, y^2]) == MonomialVector([y])
end

@testset "Random SOS should be SOS with $solver" for solver in sdp_solvers
    @polyvar x y
    x = [1, x, y, x^2, y^2, x*y]
    @test_throws ArgumentError randsos(x, monotype=:Unknown)
    for i in 1:10
        for monotype in [:Classic, :Gram]
            p = randsos(x, monotype=monotype)

            m = SOSModel(solver = solver)

            @constraint m p >= 0

            status = solve(m)

            @test status == :Optimal
        end
    end
end
