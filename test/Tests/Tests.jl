module Tests

include("utilities.jl")

include("term_test.jl")
include("bivariate_quadratic_test.jl")
include("horn_test.jl")

const linear_tests = Dict{String, Function}(
    "dsos_term" => dsos_term_test,
    "dsos_bivariate_quadratic" => dsos_bivariate_quadratic_test,
    "dsos_horn" => dsos_horn_test
)

const soc_tests = Dict{String, Function}(
    "sdsos_term" => sdsos_term_test,
    "sdsos_bivariate_quadratic" => sdsos_bivariate_quadratic_test,
    "sdsos_horn" => sdsos_horn_test
)

const sd_tests = Dict{String, Function}(
    "sos_term" => sos_term_test,
    "sos_bivariate_quadratic" => sos_bivariate_quadratic_test,
    "sos_horn" => sos_horn_test,
)

include("concave_then_convex_cubic.jl")
include("lyapunov_switched_system.jl")
include("BPT12e399.jl")

@test_suite linear
@test_suite soc
@test_suite sd

end
