@testset "CSP_graph" begin
    @testset "csp_graph" begin
        @polyvar x y z
        p = x*y + y*z
        q1 = 1 - x^2 - y^2
        q2 = 1 - y^2 - z^2
        G = csp_graph(p, [q1, q2])
        @test sort(G.nodes) == [z, y, x]
        @test length(G.edges) == 4
    end
    @testset "chordal_csp" begin
        @polyvar x y z
        p = x*y + y*z
        q1 = 1 - x^2 - y^2
        q2 = 1 - y^2 - z^2
        H, cliques = chordal_csp_graph(p, [q1, q2])
        @test length(cliques) == 2
        svar = cliques[1]∪cliques[2]
        @test CEG.contains(svar, H.nodes)
        @test CEG.contains(H.nodes, svar)
        @test sort!(svar)== sort!(H.nodes)
        @test length(H.edges) == 4
        I, cliquesI = chordal_csp_graph(p)
        @test sort!(I.nodes) == sort!(H.nodes)
        @test sort!(I.edges) ==sort!(H.edges)      
    end
end



