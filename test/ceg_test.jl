using Test

using SumOfSquares
const CEG = SumOfSquares.ChordalExtensionGraph

function fill_in_3_nodes_test(G, x, y, z)
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
    CEG.add_edge!(G, x, y)
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
    CEG.add_edge!(G, y, z)
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 1
    @test CEG.fill_in(G, z) == 0
    @testset "Disable y" begin
        H = copy(G)
        CEG.disable_node!(H, y)
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, z) == 0
        CEG.add_edge!(H, z, x)
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, z) == 0
    end
    CEG.add_edge!(G, z, x)
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
end

function fill_in_4_nodes_test(G, w, x, y, z)
    @test CEG.fill_in(G, w) == 0
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
    CEG.add_edge!(G, w, x)
    @test CEG.fill_in(G, w) == 0
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
    CEG.add_edge!(G, y, z)
    @test CEG.fill_in(G, w) == 0
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
    CEG.add_edge!(G, w, y)
    @test CEG.fill_in(G, w) == 1
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 1
    @test CEG.fill_in(G, z) == 0
    @testset "wz first" begin
        H = copy(G)
        CEG.add_edge!(H, w, z)
        @test CEG.fill_in(H, w) == 2
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, y) == 0
        @test CEG.fill_in(H, z) == 0
        @testset "Disable w" begin
            I = copy(H)
            CEG.add_edge!(I, x, z)
            @test CEG.fill_in(I, x) == 0
            @test CEG.fill_in(I, y) == 0
            @test CEG.fill_in(I, z) == 1
            CEG.add_edge!(I, x, y)
            @test CEG.fill_in(I, x) == 0
            @test CEG.fill_in(I, y) == 0
            @test CEG.fill_in(I, z) == 0
        end
        CEG.add_edge!(H, x, z)
        @test CEG.fill_in(H, w) == 1
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, y) == 0
        @test CEG.fill_in(H, z) == 1
        CEG.add_edge!(H, x, y)
        @test CEG.fill_in(H, w) == 0
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, y) == 0
        @test CEG.fill_in(H, z) == 0
    end
    CEG.add_edge!(G, x, z)
    @test CEG.fill_in(G, w) == 1
    @test CEG.fill_in(G, x) == 1
    @test CEG.fill_in(G, y) == 1
    @test CEG.fill_in(G, z) == 1
    @testset "xy then wz" begin
        H = copy(G)
        CEG.add_edge!(H, x, y)
        @test CEG.fill_in(H, w) == 0
        @test CEG.fill_in(H, x) == 1
        @test CEG.fill_in(H, y) == 1
        @test CEG.fill_in(H, z) == 0
        CEG.add_edge!(H, w, z)
        @test CEG.fill_in(H, w) == 0
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, y) == 0
        @test CEG.fill_in(H, z) == 0
    end
    CEG.add_edge!(G, w, z)
    @test CEG.fill_in(G, w) == 1
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 1
    @testset "Disable w" begin
        H = copy(G)
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, y) == 0
        @test CEG.fill_in(H, z) == 1
        CEG.add_edge!(H, x, y)
        @test CEG.fill_in(H, x) == 0
        @test CEG.fill_in(H, y) == 0
        @test CEG.fill_in(H, z) == 0
    end
    CEG.add_edge!(G, x, y)
    @test CEG.fill_in(G, w) == 0
    @test CEG.fill_in(G, x) == 0
    @test CEG.fill_in(G, y) == 0
    @test CEG.fill_in(G, z) == 0
end

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

    @testset "fill_in with 3 nodes" begin
        G = CEG.Graph()
        x = CEG.add_node!(G)
        y = CEG.add_node!(G)
        z = CEG.add_node!(G)
        fill_in_3_nodes_test(copy(G), x, y, z)
        #fill_in_3_nodes_test(CEG.FillInCache(copy(G)), x, y, z)
    end

    @testset "fill_in with 4 nodes" begin
        G = CEG.Graph()
        w = CEG.add_node!(G)
        x = CEG.add_node!(G)
        y = CEG.add_node!(G)
        z = CEG.add_node!(G)
        fill_in_4_nodes_test(copy(G), w, x, y, z)
        #fill_in_4_nodes_test(CEG.FillInCache(copy(G)), w, x, y, z)
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
