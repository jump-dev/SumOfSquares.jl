module TestDummy

using Test
import MathOptInterface as MOI
import MultivariateBases as MB
using DynamicPolynomials
using SumOfSquares

# `DummyMosek` is a stub `MOI.AbstractOptimizer` that mirrors the native
# support pattern of Mosek (see `MosekTools.jl/src/constraint.jl`):
# `MOI.PositiveSemidefiniteConeTriangle` is supported as constrained variables
# (Mosek's `barvar` matrix variables) and `MOI.Scaled{PSDConeTriangle}` is
# supported as a `VAF` constraint (Mosek's `ACC` affine conic constraints).
# Only the `MOI.supports_*` methods are implemented which is enough for the
# `LazyBridgeOptimizer` to build its cost graph and answer
# `MOI.Bridges.bridging_cost` queries.
mutable struct DummyMosek{T} <: MOI.AbstractOptimizer end

DummyMosek() = DummyMosek{Float64}()

MOI.is_empty(::DummyMosek) = true
MOI.empty!(::DummyMosek) = nothing

# Mosek `barvar` matrix variables.
function MOI.supports_add_constrained_variables(
    ::DummyMosek,
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
)
    return true
end

function MOI.supports_constraint(
    ::DummyMosek{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{
        <:Union{
            MOI.SecondOrderCone,
            MOI.RotatedSecondOrderCone,
            MOI.ExponentialCone,
            MOI.DualExponentialCone,
            MOI.PowerCone{T},
            MOI.DualPowerCone{T},
        },
    },
) where {T}
    return true
end

# Mosek `ACC` affine conic constraints.
function MOI.supports_constraint(
    ::DummyMosek{T},
    ::Type{MOI.VectorAffineFunction{T}},
    ::Type{
        <:Union{
            MOI.Zeros,
            MOI.Nonnegatives,
            MOI.SecondOrderCone,
            MOI.RotatedSecondOrderCone,
            MOI.ExponentialCone,
            MOI.DualExponentialCone,
            MOI.PowerCone{T},
            MOI.DualPowerCone{T},
            MOI.GeometricMeanCone,
            MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
        },
    },
) where {T}
    return true
end

function MOI.supports_constraint(
    ::DummyMosek{T},
    ::Type{MOI.VariableIndex},
    ::Type{
        <:Union{
            MOI.LessThan{T},
            MOI.GreaterThan{T},
            MOI.EqualTo{T},
            MOI.Interval{T},
        },
    },
) where {T}
    return true
end

function MOI.supports_constraint(
    ::DummyMosek{T},
    ::Type{MOI.ScalarAffineFunction{T}},
    ::Type{
        <:Union{
            MOI.LessThan{T},
            MOI.GreaterThan{T},
            MOI.EqualTo{T},
            MOI.Interval{T},
        },
    },
) where {T}
    return true
end

MOI.supports(::DummyMosek, ::MOI.ObjectiveSense) = true

function MOI.supports(
    ::DummyMosek{T},
    ::MOI.ObjectiveFunction{
        <:Union{MOI.VariableIndex,MOI.ScalarAffineFunction{T}},
    },
) where {T}
    return true
end

function runtests()
    for name in names(@__MODULE__; all=true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

# Verify that with `DummyMosek` (a Mosek-like stub that natively supports
# `PositiveSemidefiniteConeTriangle` as constrained variables) the bridges
# added by `SumOfSquares.Bridges.add_all_bridges` route a `WeightedSOSCone`
# constrained variable through `KernelBridge`. After removing `KernelBridge`
# the only remaining route is the `ImageBridge` constraint-side fallback,
# which is strictly more expensive.
function test_dummymosek_uses_kernel_bridge()
    T = Float64
    @polyvar x y
    optimizer = MOI.Bridges.full_bridge_optimizer(DummyMosek{T}(), T)
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
        MB.SubBasis{MB.Monomial}([y^4, x * y^3, x^2 * y^2, x^3 * y, x^4]),
        [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
        [MB.algebra_element(one(T) * x^0 * y^0)],
    )
    S = typeof(set)
    # With every bridge enabled, `KernelBridge` (variable side) wins.
    @test MOI.Bridges.bridging_cost(optimizer, S) == 6.0
    @test MOI.Bridges.is_variable_bridged(optimizer, S)
    # Removing `KernelBridge` exposes the `ImageBridge` constraint-side
    # fallback, which goes through a free variable + `ImageBridge` chain.
    MOI.Bridges.remove_bridge(
        optimizer,
        SumOfSquares.Bridges.Variable.KernelBridge{T},
    )
    @test MOI.Bridges.bridging_cost(optimizer, S) == 7.0
    @test !MOI.Bridges.is_variable_bridged(optimizer, S)
    return
end

end

TestDummy.runtests()
