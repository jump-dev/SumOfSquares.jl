include("../Tests/Tests.jl")
include("utilities.jl")

using Test, JuMP

@testset "Term" begin
    include("term.jl")
end
@testset "Bivariate quadratic" begin
    include("bivariate_quadratic.jl")
end
@testset "Concave then convex cubic" begin
    include("concave_then_convex_cubic.jl")
end
@testset "Horn" begin
    include("horn.jl")
end
@testset "Lyapunov Switched System" begin
    include("lyapunov_switched_system.jl")
end
