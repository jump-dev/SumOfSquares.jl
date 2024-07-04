@testset "CSP_graph" begin
    @testset "csp_graph" begin
        @polyvar x y z
        p = x * y + y * z
        S = @set x^2 + y^2 == 1 && y^2 + z^2 == 1
        G = SumOfSquares.Certificate.Sparsity.csp_graph(
            MB.SubBasis{MB.Monomial}(MP.monomials(p)),
            S,
        )
        @test sort(G.int2n) == [z, y, x]
    end
    @testset "chordal_csp" begin
        @polyvar x y z
        p = x * y + y * z
        S = @set x^2 + y^2 == 1 && y^2 + z^2 == 1
        H, cliques = SumOfSquares.Certificate.Sparsity.chordal_csp_graph(
            MB.SubBasis{MB.Monomial}(MP.monomials(p)),
            S,
        )
        @test length(cliques) == 2
        svar = cliques[1] ∪ cliques[2]
        @test H.int2n ⊆ svar
        @test svar ⊆ H.int2n
        @test sort!(svar) == sort!(H.int2n)
        I, cliquesI = SumOfSquares.Certificate.Sparsity.chordal_csp_graph(
            MB.SubBasis{MB.Monomial}(MP.monomials(p)),
            FullSpace(),
        )
        @test sort!(I.int2n) == sort!(H.int2n)
        @test sort!(collect.(I.graph.neighbors)) ==
              sort!(collect.(H.graph.neighbors))
    end
end
