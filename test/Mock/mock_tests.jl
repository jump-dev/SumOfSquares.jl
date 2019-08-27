include("../Tests/Tests.jl")
include("utilities.jl")

using Test, JuMP

@testset "Term" begin
    include("term.jl")
end
@testset "Term fixed" begin
    include("term_fixed.jl")
end
@testset "Quadratic" begin
    include("quadratic.jl")
end
@testset "Choi" begin
    include("choi.jl")
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
@testset "Motzkin" begin
    include("motzkin.jl")
end
@testset "BPT12e399" begin
    include("BPT12e399.jl")
end
@testset "Max Cut" begin
    include("maxcut.jl")
end
@testset "Chebyshev" begin
    include("chebyshev.jl")
end
@testset "Quartic comparison" begin
    include("quartic_comparison.jl")
end
@testset "SOSDEMO9" begin
    include("sosdemo9.jl")
end
@testset "SOSDEMO10" begin
    include("sosdemo10.jl")
end
