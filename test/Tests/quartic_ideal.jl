using Test
import MultivariatePolynomials as MP
using SumOfSquares
using DynamicPolynomials

function quartic_ideal_test(
    optimizer,
    config::MOI.Test.Config,
    degree::Union{Nothing,Int},
    remainder::Bool,
)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @polyvar x
    # Set {-1, 0, 1}
    K = @set x^3 == x
    p = (x^2 - 1)^2
    cref = @constraint(
        model,
        p in SOSCone(),
        domain = K,
        maxdegree = degree,
        newton_of_remainder = remainder
    )
    optimize!(model)

    if remainder && (degree === nothing || degree < 4)
        @test termination_status(model) in
              [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    else
        @test termination_status(model) == MOI.OPTIMAL
        @test primal_status(model) == MOI.FEASIBLE_POINT
        μ = dual(cref)
        @test MP.monomial(moments(μ)[1].polynomial) == 1
        @test MP.monomial(moments(μ)[2].polynomial) == x^2
        @test MP.monomial(moments(μ)[3].polynomial) == x^4
        μ = moments(cref)
        @test MP.monomial(moments(μ)[1].polynomial) == 1
        @test MP.monomial(moments(μ)[2].polynomial) == x^1
        @test MP.monomial(moments(μ)[3].polynomial) == x^2
        @test MB.keys_as_monomials(moment_matrix(cref).basis) == [1, x, x^2]
    end
end
function quartic_ideal_test(optimizer, config)
    return quartic_ideal_test(optimizer, config, nothing, false)
end
sd_tests["quartic_ideal"] = quartic_ideal_test
function quartic_ideal_4_test(optimizer, config)
    return quartic_ideal_test(optimizer, config, 4, false)
end
sd_tests["quartic_ideal_4"] = quartic_ideal_4_test
function quartic_ideal_rem_test(optimizer, config)
    return quartic_ideal_test(optimizer, config, nothing, true)
end
sd_tests["quartic_ideal_rem"] = quartic_ideal_rem_test
function quartic_ideal_2_rem_test(optimizer, config)
    return quartic_ideal_test(optimizer, config, 2, true)
end
sd_tests["quartic_ideal_2_rem"] = quartic_ideal_2_rem_test
function quartic_ideal_4_rem_test(optimizer, config)
    return quartic_ideal_test(optimizer, config, 4, true)
end
sd_tests["quartic_ideal_4_rem"] = quartic_ideal_4_rem_test
