using Test
using SumOfSquares
using DynamicPolynomials

function quartic_ideal_test(optimizer, config::MOIT.TestConfig,
                            degree::Union{Nothing, Int}, remainder::Bool)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @polyvar x
    # Set {-1, 0, 1}
    K = @set x^3 == x
    p = (x^2 - 1)^2
    cref = @constraint(model, p in SOSCone(), domain = K, maxdegree=degree, remainder = remainder)
    optimize!(model)

    if remainder && (degree === nothing || degree < 4)
        @test termination_status(model) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    else
        @test termination_status(model) == MOI.OPTIMAL
        @test primal_status(model) == MOI.FEASIBLE_POINT
        @show dual(cref)
        @show moment_matrix(cref)
        @show moments(cref)
    end
end
quartic_ideal_test(optimizer, config)   = quartic_ideal_test(optimizer, config, nothing, false)
sd_tests["quartic_ideal"] = quartic_ideal_test
quartic_ideal_4_test(optimizer, config) = quartic_ideal_test(optimizer, config, 4, false)
sd_tests["quartic_ideal_4"] = quartic_ideal_4_test
quartic_ideal_rem_test(optimizer, config)   = quartic_ideal_test(optimizer, config, nothing, true)
sd_tests["quartic_ideal_rem"] = quartic_ideal_rem_test
quartic_ideal_2_rem_test(optimizer, config) = quartic_ideal_test(optimizer, config, 2, true)
sd_tests["quartic_ideal_2_rem"] = quartic_ideal_2_rem_test
quartic_ideal_4_rem_test(optimizer, config) = quartic_ideal_test(optimizer, config, 4, true)
sd_tests["quartic_ideal_4_rem"] = quartic_ideal_4_rem_test
