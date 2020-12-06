using Test

include("all_examples.jl")

@testset "run_examples.jl" begin
    @testset "$(example)" for example in EXAMPLES
        include(example)
    end
end
