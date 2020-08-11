using Test

const EXAMPLES = filter(ex -> endswith(ex, ".jl") && ex != "run_examples.jl" && ex != "goldsteinprice.jl" && ex != "sparse_polynomials.jl",
                        readdir(@__DIR__))

@testset "run_examples.jl" begin
    @testset "$(example)" for example in EXAMPLES
        include(example)
    end
end
