using Test

using SumOfSquares
const CEG = SumOfSquares.ChordalExtensionGraph

@testset "Chordal Extensions" begin
    @testset "CEG.Graph" begin
        @test CEG.Graph{Int}() isa CEG.Graph{Int}
        @test_throws ErrorException CEG.Graph()
        G = CEG.Graph{Int}()
        @test CEG.add_node!(G, 1) == 1
        @test CEG.add_edge!(G, (1, 2)) == (1, 2)
        @test CEG.add_edge!.(G, [(1, 2), (2, 3)]) isa Vector{Tuple{Int, Int}}
        @test CEG.add_clique!(G, [1, 2, 3]) == nothing
    end

    @testset "show" begin
        G = CEG.Graph{Symbol}()
        @test sprint(show, G) == "Symbol-Graph with\nNodes:\nEdges:\n"
        CEG.add_edge!(G, :x, :y)
        @test sprint(show, G) isa String
    end

    @testset "sub_sets" begin
        G = CEG.Graph{Symbol}()
        CEG.add_clique!(G, [:x, :y, :z])
        H = CEG.sub_graph(G, [:x, :y])
        @test H.int2n == [:x, :y]
        @test H.edges == [[2], [1]]
        @test CEG.neighbors(G, :x) == [:y, :z]
    end

    @testset "cliques" begin
        G = CEG.Graph{Symbol}()
        CEG.add_edge!(G, :x, :y)
        CEG.add_edge!(G, :y, :z)
        @test CEG.fill_in.(G, [:x, :y, :z]) == [0; 1; 0]
        @test !CEG.is_clique(G, [:x, :y, :z])
        @test CEG.is_clique(G, [:x, :y])
    end

    @testset "chordal" begin
        G = CEG.Graph{Symbol}()
        CEG.add_edge!.(G, [(:x, :y), (:y, :z)])
        H, cliques = CEG.chordal_extension(G)
        @test H.int2n == G.int2n
        @test H.edges == G.edges
        @test cliques == [[:y, :z], [:x, :y]]

        G = CEG.Graph{Int}()
        CEG.add_edge!.(G, [(1, 2), (1, 3), (3, 4), (2, 4)])
        H, cliques = CEG.chordal_extension(G)
        @test CEG.n_edges(G) <= CEG.n_edges(H)
        @test length.(cliques) == [3, 3]

    end

end
