module Tests

include("utilities.jl")

include("term_test.jl")
include("bivariate_quadratic_test.jl")
include("concave_then_convex_cubic_test.jl")
include("horn_test.jl")

const linear_tests = Dict(
    "dsos_term" => dsos_term_test,
    "dsos_bivariate_quadratic" => dsos_bivariate_quadratic_test,
    "dsos_concave_then_convex_cubic" => dsos_concave_then_convex_cubic_test,
    "dsos_horn" => dsos_horn_test
)

const soc_tests = Dict(
    "sdsos_term" => sdsos_term_test,
    "sdsos_bivariate_quadratic" => sdsos_bivariate_quadratic_test,
    "sdsos_concave_then_convex_cubic" => sdsos_concave_then_convex_cubic_test,
    "sdsos_horn" => sdsos_horn_test
)

const sd_tests = Dict(
    "sos_term" => sos_term_test,
    "sos_bivariate_quadratic" => sos_bivariate_quadratic_test,
    "sos_concave_then_convex_cubic" => sos_concave_then_convex_cubic_test,
    "sos_horn" => sos_horn_test,
)

include("lyapunov_switched_system.jl")

@test_suite linear
@test_suite soc
@test_suite sd

end
