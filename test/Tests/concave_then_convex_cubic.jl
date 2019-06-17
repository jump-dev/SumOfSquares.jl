using Test
using SumOfSquares
using DynamicPolynomials

# Inspired from https://github.com/JuliaOpt/SumOfSquares.jl/issues/79
function concave_then_convex_cubic_test(optimizer, config::MOIT.TestConfig,
                                        MCT::Type{<:MOI.AbstractVectorSet})
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @polyvar x
    @variable(model, p, Poly(monomials(x, 0:3)))
    cone = ConvexPolyInnerCone{MCT}()
    cref_convex  = @constraint(model,  p in cone, domain = (@set x >= 0))
    cref_concave = @constraint(model, -p in cone, domain = (@set x <= 0))
    @constraint(model, p(x =>  2) ==  8)
    @constraint(model, p(x => -2) == -8)
    @constraint(model, p(x =>  1) ==  1)
    @constraint(model, p(x => -1) == -1)

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(p) ≈ x^3 atol=atol rtol=rtol
    @test_throws SumOfSquares.ValueNotSupported value(cref_convex)
    @test_throws SumOfSquares.ValueNotSupported value(cref_concave)

    # The monomials contain the variables created for the Hessian so we cannot
    # check them easily
    for μ in [dual(cref_convex), dual(cref_concave)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 2
    end
    for μ in [moments(cref_convex), moments(cref_concave)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 7
    end
end

function sos_concave_then_convex_cubic_test(optimizer, config)
    concave_then_convex_cubic_test(optimizer, config,
                                   MOI.PositiveSemidefiniteConeTriangle)
end
sd_tests["sos_concave_then_convex_cubic"] = sos_concave_then_convex_cubic_test
function sdsos_concave_then_convex_cubic_test(optimizer, config)
    concave_then_convex_cubic_test(optimizer, config,
                                   SumOfSquares.ScaledDiagonallyDominantConeTriangle)
end
soc_tests["sdsos_concave_then_convex_cubic"] = sdsos_concave_then_convex_cubic_test
function dsos_concave_then_convex_cubic_test(optimizer, config)
    concave_then_convex_cubic_test(optimizer, config,
                                   SumOfSquares.DiagonallyDominantConeTriangle)
end
linear_tests["dsos_concave_then_convex_cubic"] = dsos_concave_then_convex_cubic_test
