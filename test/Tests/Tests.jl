module Tests

include("utilities.jl")

const linear_tests = Dict{String, Function}()
const soc_tests = Dict{String, Function}()
const sd_tests = Dict{String, Function}()

include("term.jl")
include("term_fixed.jl")
include("quadratic.jl")
include("choi.jl")
include("horn.jl")
include("concave_then_convex_cubic.jl")
include("lyapunov_switched_system.jl")
include("motzkin.jl")
include("BPT12e399.jl")
include("maxcut.jl")
include("chebyshev.jl")
include("quartic_comparison.jl")
include("sosdemo9.jl")
include("sosdemo10.jl")

@test_suite linear
@test_suite soc
@test_suite sd

end
