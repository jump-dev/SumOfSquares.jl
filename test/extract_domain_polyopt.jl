# Test for the example https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb

using JuMP
using SumOfSquares
using SemialgebraicSets
using MultivariateMoments

@testset "Polynomial Optimization example with $(typeof(solver))" for solver in sdp_solvers
    isscs(solver) && continue
    @polyvar x y
    p = x^3 - x^2 + 2x*y -y^2 + y^3
    S = @set x >= 0 && y >= 0 && x + y >= 1

    for (maxdeg, found) in [(3, false), (4, true), (5, true)]
        MOI.empty!(solver)
        m = SOSModel(optimizer=solver)
        @variable m α
        @objective m Max α
        c = @constraint m p >= α domain = S maxdegree = maxdeg

        JuMP.optimize(m)
        @test JuMP.terminationstatus(m) == MOI.Success
        @test JuMP.objectivevalue(m) ≈ 0 atol=1e-4

        μ = JuMP.resultdual(c)
        X = certificate_monomials(PolyJuMP.getdelegate(c))
        ν = matmeasure(μ, X)
        ranktol = 1e-3
        atoms = extractatoms(ν, ranktol)
        @test isnull(atoms) == !found
        if !isnull(atoms)
            η = get(atoms)
            @test η.atoms[1].weight ≈ 1/2 atol=1e-2
            @test η.atoms[2].weight ≈ 1/2 atol=1e-2
            @test isapprox(η.atoms[1].center, [0, 1], atol=1e-2) || isapprox(η.atoms[1].center, [1, 0], atol=1e-2)
            @test isapprox(η.atoms[2].center, [0, 1], atol=1e-2) || isapprox(η.atoms[2].center, [1, 0], atol=1e-2)
        end
    end
end
