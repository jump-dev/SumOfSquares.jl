@testset "CSP_graph" begin
    @testset "csp_graph" begin
        @polyvar x y z
        p = x*y + y*z
        q1 = 1 - x^2 - y^2
        q2 = 1 - y^2 - z^2
        G = csp_graph(p, [q1, q2])
        @test sort(G.int2n) == [z, y, x]
        @test CEG.n_edges(G) == 4
    end
    @testset "chordal_csp" begin
        @polyvar x y z
        p = x*y + y*z
        q1 = 1 - x^2 - y^2
        q2 = 1 - y^2 - z^2
        H, cliques = chordal_csp_graph(p, [q1, q2])
        @test length(cliques) == 2
        svar = cliques[1]∪cliques[2]
        @test H.int2n ⊆ svar
        @test svar ⊆ H.int2n
        @test sort!(svar)== sort!(H.int2n)
        @test CEG.n_edges(H) == 4
        I, cliquesI = chordal_csp_graph(p)
        @test sort!(I.int2n) == sort!(H.int2n)
        @test sort!(I.edges) == sort!(H.edges)
    end
end
