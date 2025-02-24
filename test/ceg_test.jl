using Test

using SumOfSquares
const CEG = SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph


@testset "Chordal Extensions" begin
    @testset "CEG.LabelledGraph" begin
        @test CEG.LabelledGraph{Int}() isa CEG.LabelledGraph{Int}
        @test_throws ErrorException CEG.LabelledGraph()
        G = CEG.LabelledGraph{Int}()
        @test CEG.add_node!(G, 1) == 1
        @test CEG.add_edge!(G, (1, 2)) == (1, 2)
        @test CEG.add_edge!.(G, [(1, 2), (2, 3)]) isa Vector{Tuple{Int,Int}}
        @test CEG.add_clique!(G, [1, 2, 3]) === nothing
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
        @test !CEG.is_clique(G.graph, [1, 2, 3])
        @test CEG.is_clique(G.graph, [1, 2])
    end

    @testset "chordal" begin
        G = CEG.LabelledGraph{Symbol}()
        CEG.add_edge!.(G, [(:x, :y), (:y, :z)])
        H, cliques = CEG.chordal_extension(G, CEG.GreedyFillIn())
        @test H.int2n == G.int2n
        @test H.graph.neighbors == G.graph.neighbors
        @test cliques == [[:z, :y], [:x, :y]]

        G = CEG.LabelledGraph{Int}()
        CEG.add_edge!.(G, [(1, 2), (1, 3), (3, 4), (2, 4)])
        H, cliques = CEG.chordal_extension(G, CEG.GreedyFillIn())
        @test length.(cliques) == [3, 3]
    end

    @testset "Add edges between two `chordal_extension` calls" begin
        G = CEG.LabelledGraph{Int}()
        CEG.add_node!.(G, [1, 2, 3, 4, 5])
        CEG.add_edge!(G, 1, 2)
        CEG.add_edge!(G, 2, 3)
        CEG.add_edge!(G, 3, 4)
        CEG.add_edge!(G, 1, 4)

        H, cliques = CEG.chordal_extension(G, CEG.GreedyFillIn())
        @test cliques == [[4, 1, 3], [2, 1, 3], [5]]

        CEG.add_edge!(H, 1, 5)
        I, cliques = CEG.chordal_extension(H, CEG.GreedyFillIn())
        @test cliques == [[5, 1], [4, 1, 3], [2, 1, 3]]
    end
end
