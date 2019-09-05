@testset "Chordal Extensions" begin
    @testset "CEG.Graph" begin
        @test CEG.Graph{Int}() isa CEG.Graph{Int}
        @test CEG.Graph() == nothing
        G = CEG.Graph{Int}()
        @test CEG.add_node!(G, 1) == 1
        @test CEG.add_edge!(G, (1, 2)) == Vector{Tuple{Int, Int}}([(1, 2), (2, 1)])
        @test CEG.add_edge!.(G, [(1, 2), (2, 3)]) isa Vector{Vector{Tuple{Int, Int}}}
        @test CEG.add_clique!(G, [1, 2, 3]) == nothing
    end

    @testset "show" begin
        G = CEG.Graph{Symbol}()
        @test sprint(show, G) == "Symbol-Graph with\nNodes:\nEdges:\n"
        CEG.add_edge!(G, :x, :y)
        @test sprint(show, G) == "Symbol-Graph with\nNodes:\nx, y, Edges:\n(:x, :y)\n(:y, :x)\n"
    end

    @testset "sub_sets" begin
        G = CEG.Graph{Symbol}()
        CEG.add_clique!(G, [:x, :y, :z])
        H = CEG.sub_graph(G, [:x, :y])
        @test H.nodes == [:x, :y]
        @test H.edges == [(:x, :y), (:y, :x)]
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

    @testset "contains" begin
        foo  = [1, 2, 3]
        bar = [1, 2]
        goo = [1, 2, 4]
        mar = [3, 2]
        @test CEG.contains(foo, foo)
        @test CEG.contains(foo, bar)
        @test !CEG.contains(foo, goo)    
        @test CEG.contains(foo, mar)
        @test !CEG.contains(bar, foo)
        @test !CEG.contains(bar, goo)
        @test !CEG.contains(bar, mar)
        @test !CEG.contains(goo, mar)
    end


    @testset "chordal" begin
        G = CEG.Graph{Symbol}()
        CEG.add_edge!.(G, [(:x, :y), (:y, :z)])
        H, cliques = CEG.chordal_extension(G)
        @test H.nodes == G.nodes
        @test H.edges == G.edges
        @test cliques == [[:y, :z], [:x, :y]]

        G = CEG.Graph{Int}()
        CEG.add_edge!.(G, [(1, 2), (1, 3), (3, 4), (2, 4)])
        H, cliques = CEG.chordal_extension(G)
        @test length(G.edges) < length(H.edges)
        @test length.(cliques) == [3, 3]
    end

    @testset "adjacency_matrix" begin
        G = CEG.Graph{Int}()
        CEG.add_edge!.(G, [(1, 2), (2, 3)])
        @test CEG.adjacency_matrix(G) isa SparseArrays.SparseMatrixCSC{Int,Int}
        A = CEG.adjacency_matrix(G)
        @test nnz(A) == 7
        @test sprint(show, Matrix(A)) == "[1 1 0; 1 1 1; 0 1 1]"
    end
end
