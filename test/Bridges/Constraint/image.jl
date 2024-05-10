module TestConstraintImage

using Test
using DynamicPolynomials
using SumOfSquares

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

# [Parrilo2003, Example 6.1]
function test_runtests()
    T = Float64
    @polyvar x y
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.ImageBridge,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(
                    MOI.VectorAffineTerm{T}[],
                    T[5, -1, 2, 2],
                ),
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MonomialBasis([y^4, x^2 * y^2, x^3 * y, x^4]),
                    [MonomialBasis([y^2, x * y, x^2])],
                    [one(T) * x^0 * y^0],
                ),
            )
        end,
        model -> begin
            λ = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.Utilities.operate(
                    vcat,
                    T,
                    T(5),
                    T(0),
                    T(2) * λ - T(1),
                    T(-1) * λ,
                    T(1),
                    T(2),
                ),
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
        end,
    )
    return
end

end  # module

TestConstraintImage.runtests()
