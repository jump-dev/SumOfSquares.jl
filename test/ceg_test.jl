using Test

using SumOfSquares
const CEG = SumOfSquares.ChordalExtensionGraph

@testset "Chordal Extensions" begin
    @testset "CEG.LabelledGraph" begin
        @test CEG.LabelledGraph{Int}() isa CEG.LabelledGraph{Int}
        @test_throws ErrorException CEG.LabelledGraph()
        G = CEG.LabelledGraph{Int}()
        @test CEG.add_node!(G, 1) == 1
        @test CEG.add_edge!(G, (1, 2)) == (1, 2)
        @test CEG.add_edge!.(G, [(1, 2), (2, 3)]) isa Vector{Tuple{Int, Int}}
        @test CEG.add_clique!(G, [1, 2, 3]) == nothing
    end

    @testset "show" begin
        G = CEG.LabelledGraph{Symbol}()
        @test sprint(show, G) == "Symbol-LabelledGraph with\nNodes:\nEdges:\n"
        CEG.add_edge!(G, :x, :y)
        @test sprint(show, G) isa String
    end

    @testset "sub_sets" begin
        G = CEG.LabelledGraph{Symbol}()
        CEG.add_clique!(G, [:x, :y, :z])
        @test CEG.neighbors(G.graph, 1) == Set([2, 3])
    end

    @testset "cliques" begin
        G = CEG.LabelledGraph{Symbol}()
        CEG.add_edge!(G, :x, :y)
        CEG.add_edge!(G, :y, :z)
        @test CEG.fill_in.(G.graph, [1, 2, 3]) == [0, 1, 0]
        @test !CEG.is_clique(G.graph, [1, 2, 3])
        @test CEG.is_clique(G.graph, [1, 2])
    end

    @testset "chordal" begin
        G = CEG.LabelledGraph{Symbol}()
        CEG.add_edge!.(G, [(:x, :y), (:y, :z)])
        H, cliques = CEG.chordal_extension(G)
        @test H.int2n == G.int2n
        @test H.graph.neighbors == G.graph.neighbors
        @test cliques == [[:y, :z], [:x, :y]]

        G = CEG.LabelledGraph{Int}()
        CEG.add_edge!.(G, [(1, 2), (1, 3), (3, 4), (2, 4)])
        H, cliques = CEG.chordal_extension(G)
        @test length.(cliques) == [3, 3]

    end

end
