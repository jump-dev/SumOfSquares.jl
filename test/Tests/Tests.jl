module Tests

include("utilities.jl")

const linear_tests = Dict{String,Function}()
const soc_tests = Dict{String,Function}()
const sd_tests = Dict{String,Function}()

include("term.jl")
include("term_fixed.jl")
include("quartic_constant.jl")
include("quadratic.jl")
include("quartic_ideal.jl")
include("univariate_sum.jl")
include("rearrangement.jl")
include("choi.jl")
include("simple_matrix.jl")
include("horn.jl")
include("concave_then_convex_cubic.jl")
include("lyapunov_switched_system.jl")
include("motzkin.jl")
include("BPT12e399.jl")
include("maxcut.jl")
include("chebyshev.jl")
include("quartic_comparison.jl")
include("options_pricing.jl")
include("sosdemo5.jl")
include("sosdemo9.jl")
include("sosdemo10.jl")

@test_suite linear
@test_suite soc
@test_suite sd

end
