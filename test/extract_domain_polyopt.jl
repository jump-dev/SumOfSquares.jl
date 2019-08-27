# Test for the example https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb

using JuMP
using SumOfSquares
using MultivariateMoments

@testset "Polynomial Optimization example with $(factory.constructor)" for factory in sdp_factories
    isscs(factory) && continue
    @polyvar x y
    p = x^3 - x^2 + 2x*y -y^2 + y^3
    S = @set x >= 0 && y >= 0 && x + y >= 1

    for (maxdeg, found) in [(3, false), (4, true), (5, true)]
        m = SOSModel(factory)
        @variable m α
        # TODO Change 1α to α once we have variable bridges
        @objective m Max 1α
        c = @constraint m p >= α domain = S maxdegree = maxdeg

        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.objective_value(m) ≈ 0 atol=1e-4

        for λ in lagrangian_multipliers(c)
            @test all(eigvals(Matrix(λ.Q)) .>= -1e-2)
        end

        X = certificate_monomials(c)
        ν = moment_matrix(c)
        ranktol = 1e-3
        atoms = extractatoms(ν, ranktol)
        @test (atoms === nothing) == !found
        if atoms !== nothing
            η = atoms
            @test η.atoms[1].weight ≈ 1/2 atol=1e-2
            @test η.atoms[2].weight ≈ 1/2 atol=1e-2
            @test isapprox(η.atoms[1].center, [0, 1], atol=1e-2) || isapprox(η.atoms[1].center, [1, 0], atol=1e-2)
            @test isapprox(η.atoms[2].center, [0, 1], atol=1e-2) || isapprox(η.atoms[2].center, [1, 0], atol=1e-2)
        end
    end
end
