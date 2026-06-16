# Copyright (c) 2026: Benoît Legat and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.
#
# Self-contained tests for the `Bridges.Variable.LowRankBridge` parameterisation
# rewrite (the change that lets `MultivariateBases.TrigEvalMatrix` survive
# intact to `LowRankOpt.BurerMonteiro`).
#
# What this file pins down:
#
# 1. `LowRankBridge` now carries the basis / gram-basis / weight type
#    parameters `{T,M,B,G,W}` (was `{T,M}`).
# 2. `_transformation_type(G, B)` / `_row_view_type(MT)` resolve at the
#    type level (via `Base.promote_op`) so `added_constrained_variable_types`
#    declares the concrete SubArray factor type that
#    `bridge_constrained_variable` actually produces.
# 3. The bridge stores `view(U, j, :)` (not `collect(...)`) inside each
#    `LRO.Factorization`, so `parent` on the factor returns the underlying
#    matrix `U` — the contract that the (future) batched-FFT path in
#    BurerMonteiro relies on.
# 4. For `gram_basis = MultivariateBases.Trigonometric`, the same path
#    resolves `U::MultivariateBases.TrigEvalMatrix` and the factor type is
#    a `SubArray` of that matrix.

module TestVariableLowRankParam

using Test
using DynamicPolynomials
using JuMP
using SumOfSquares
import MultivariateBases as MB
import LowRankOpt as LRO
import StarAlgebras as SA
import MathOptInterface as MOI

# Build a `WeightedSOSCone{M,B,G,W}` of the shape `LowRankBridge` consumes,
# parameterised by the gram-basis type so we can exercise both the dense
# (`Monomial`) and FFT (`Trigonometric`) paths.
function _weighted_sos_cone(::Type{B}; var = :t, gram_degree = 2, n_pts = 5) where {B}
    @polyvar t
    pts = [Float64[k / (n_pts - 1)] for k in 0:(n_pts-1)]
    lag = MB.LagrangeBasis((t,), pts)
    gram = MB.SubBasis{B}(monomials(t, 0:gram_degree))
    weight = MB.algebra_element(1.0 * t^0)
    return SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
        lag,
        [gram],
        [weight],
    ), lag, gram, weight
end

# 1 ─────────────────────────────────────────────────────────────────────────
# Bridge type carries 5 parameters.
# ─────────────────────────────────────────────────────────────────────────

function test_bridge_struct_has_five_type_parameters()
    BridgeT = SumOfSquares.Bridges.Variable.LowRankBridge
    @test BridgeT isa UnionAll
    # `LowRankBridge{T}` should still be a UnionAll over the remaining 4
    # parameters — the existing bridge-graph registration uses that syntax.
    @test BridgeT{Float64} isa UnionAll
    # Concretely instantiate to make sure 5 type parameters resolve.
    set, _, gram, weight = _weighted_sos_cone(MB.Monomial)
    SetT = typeof(set)
    Mp = MOI.PositiveSemidefiniteConeTriangle
    Bp = typeof(set.basis)
    Gp = eltype(set.gram_bases)
    Wp = eltype(set.weights)
    concrete = BridgeT{Float64,Mp,Bp,Gp,Wp}
    @test concrete <: BridgeT
    @test isconcretetype(concrete) || concrete isa DataType
end

# 2 ─────────────────────────────────────────────────────────────────────────
# `concrete_bridge_type` round-trips: given a `WeightedSOSCone{M,B,G,W}`
# it returns the 5-param `LowRankBridge{T,M,B,G,W}`.
# ─────────────────────────────────────────────────────────────────────────

function test_concrete_bridge_type_monomial_gram()
    set, _, _, _ = _weighted_sos_cone(MB.Monomial)
    bridge_t = MOI.Bridges.Variable.concrete_bridge_type(
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
        typeof(set),
    )
    @test bridge_t isa Type
    @test bridge_t <: SumOfSquares.Bridges.Variable.LowRankBridge{Float64}
    @test length(bridge_t.parameters) == 5
    @test bridge_t.parameters[1] === Float64
    @test bridge_t.parameters[2] === MOI.PositiveSemidefiniteConeTriangle
    @test bridge_t.parameters[3] === typeof(set.basis)
    @test bridge_t.parameters[4] === eltype(set.gram_bases)
    @test bridge_t.parameters[5] === eltype(set.weights)
end

function test_concrete_bridge_type_trigonometric_gram()
    set, _, _, _ = _weighted_sos_cone(MB.Trigonometric; gram_degree = 4, n_pts = 9)
    bridge_t = MOI.Bridges.Variable.concrete_bridge_type(
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
        typeof(set),
    )
    @test bridge_t <: SumOfSquares.Bridges.Variable.LowRankBridge{Float64}
    @test bridge_t.parameters[4] === eltype(set.gram_bases)
end

# 3 ─────────────────────────────────────────────────────────────────────────
# `_transformation_type(G, B)` is a compile-time map from
# `(gram_basis_type, target_basis_type)` to the matrix type that
# `MB.transformation_to` produces. This is what
# `added_constrained_variable_types` consults.
# ─────────────────────────────────────────────────────────────────────────

function test_transformation_type_monomial_returns_matrix()
    _, lag, gram, _ = _weighted_sos_cone(MB.Monomial)
    MT = SumOfSquares.Bridges.Variable._transformation_type(typeof(gram), typeof(lag))
    @test MT !== Any
    @test MT <: AbstractMatrix{Float64}
    @test MT === typeof(MB.transformation_to(gram, lag))
end

function test_transformation_type_trigonometric_returns_trig_eval()
    _, lag, gram, _ = _weighted_sos_cone(MB.Trigonometric; gram_degree = 4, n_pts = 9)
    MT = SumOfSquares.Bridges.Variable._transformation_type(typeof(gram), typeof(lag))
    @test MT <: MB.TrigEvalMatrix{Float64}
    @test MT === typeof(MB.transformation_to(gram, lag))
end

function test_row_view_type_matches_view_at_runtime()
    _, lag, gram, _ = _weighted_sos_cone(MB.Monomial)
    MT = SumOfSquares.Bridges.Variable._transformation_type(typeof(gram), typeof(lag))
    FT = SumOfSquares.Bridges.Variable._row_view_type(MT)
    U = MB.transformation_to(gram, lag)
    @test FT === typeof(view(U, 1, :))
end

function test_row_view_type_for_trig_eval()
    _, lag, gram, _ = _weighted_sos_cone(MB.Trigonometric; gram_degree = 4, n_pts = 9)
    MT = SumOfSquares.Bridges.Variable._transformation_type(typeof(gram), typeof(lag))
    FT = SumOfSquares.Bridges.Variable._row_view_type(MT)
    U = MB.transformation_to(gram, lag)
    @test FT === typeof(view(U, 1, :))
    # The factor is a `SubArray` whose parent IS the `TrigEvalMatrix` — the
    # invariant downstream batched-FFT consumers will rely on.
    v = view(U, 2, :)
    @test parent(v) === U
    @test parent(v) isa MB.TrigEvalMatrix
end

# 4 ─────────────────────────────────────────────────────────────────────────
# `added_constrained_variable_types` declares the *same* concrete
# `SetDotProducts` type that `bridge_constrained_variable` actually adds.
# (This is the contract MOI's bridge graph relies on; the older `Vector{T}`
# declaration here would silently mismatch the new `SubArray` factor type
# and PolyJuMP would refuse to bridge.)
# ─────────────────────────────────────────────────────────────────────────

function _added_constraint_set_type(BridgeT)
    types = MOI.Bridges.added_constrained_variable_types(BridgeT)
    @assert length(types) == 1
    return types[1][1]
end

function test_added_constrained_variable_types_monomial()
    set, _, _, _ = _weighted_sos_cone(MB.Monomial)
    BridgeT = MOI.Bridges.Variable.concrete_bridge_type(
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
        typeof(set),
    )
    declared = _added_constraint_set_type(BridgeT)
    @test declared <: LRO.SetDotProducts{LRO.WITHOUT_SET}
    # 4 type parameters (W, S, V, Vs) — the relaxed `SetDotProducts`.
    @test length(declared.parameters) == 4
    # Factor type matches `view(::Matrix{Float64}, ::Int, ::Colon)`.
    V = declared.parameters[3]
    @test V <: LRO.TriangleVectorization{Float64}
    Fact = V.parameters[2]
    @test Fact <: LRO.Factorization{Float64}
    F = Fact.parameters[2]
    @test F === typeof(view(Matrix{Float64}(undef, 1, 1), 1, :))
end

function test_added_constrained_variable_types_trigonometric()
    set, lag, gram, _ = _weighted_sos_cone(
        MB.Trigonometric;
        gram_degree = 4,
        n_pts = 9,
    )
    BridgeT = MOI.Bridges.Variable.concrete_bridge_type(
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
        typeof(set),
    )
    declared = _added_constraint_set_type(BridgeT)
    @test declared <: LRO.SetDotProducts{LRO.WITHOUT_SET}
    V = declared.parameters[3]
    Fact = V.parameters[2]
    F = Fact.parameters[2]
    # Factor type must be a `SubArray` *of* the `TrigEvalMatrix`, not of a
    # plain `Matrix{Float64}` — this is what preserves the FFT path through
    # the bridge chain.
    U = MB.transformation_to(gram, lag)
    @test F === typeof(view(U, 1, :))
    @test F <: SubArray
    @test parent(view(U, 1, :)) isa MB.TrigEvalMatrix
end

# 5 ─────────────────────────────────────────────────────────────────────────
# End-to-end: after `bridge_constrained_variable`, the constraint that lands
# at the inner model has the declared concrete type, and pulling out a
# constraint vector gives a `TriangleVectorization` whose `Factorization`'s
# `.factor` is a `SubArray` with a recoverable parent matrix.
# ─────────────────────────────────────────────────────────────────────────

function _bridge_into_inner(::Type{B}; kwargs...) where {B}
    set, _, _, _ = _weighted_sos_cone(B; kwargs...)
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Variable.SingleBridgeOptimizer{
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
    }(inner)
    MOI.add_constrained_variables(model, set)
    return inner, set
end

function test_bridge_emits_declared_set_type_monomial()
    inner, set = _bridge_into_inner(MB.Monomial)
    BridgeT = MOI.Bridges.Variable.concrete_bridge_type(
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
        typeof(set),
    )
    declared = _added_constraint_set_type(BridgeT)
    cts = MOI.get(inner, MOI.ListOfConstraintTypesPresent())
    @test any(((F, S),) -> S === declared, cts)
end

function test_bridge_emits_declared_set_type_trigonometric()
    inner, set = _bridge_into_inner(MB.Trigonometric; gram_degree = 4, n_pts = 9)
    BridgeT = MOI.Bridges.Variable.concrete_bridge_type(
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
        typeof(set),
    )
    declared = _added_constraint_set_type(BridgeT)
    cts = MOI.get(inner, MOI.ListOfConstraintTypesPresent())
    @test any(((F, S),) -> S === declared, cts)
end

function test_bridge_stores_subarray_factor_with_recoverable_parent()
    # The bridge must use `view(U, j, :)` rather than materialising the row
    # so that `parent(.factor) === U` for each `j`. This is what lets the
    # downstream batched-FFT path detect "all rank-1 constraints share a
    # common parent matrix".
    inner, set = _bridge_into_inner(MB.Trigonometric; gram_degree = 4, n_pts = 9)
    declared = _added_constraint_set_type(
        MOI.Bridges.Variable.concrete_bridge_type(
            SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
            typeof(set),
        ),
    )
    cis = MOI.get(
        inner,
        MOI.ListOfConstraintIndices{MOI.VectorOfVariables,declared}(),
    )
    @test length(cis) == 1
    s = MOI.get(inner, MOI.ConstraintSet(), first(cis))
    @test length(s.vectors) == length(set.basis)
    factors = [tv.matrix.factor for tv in s.vectors]
    @test all(f -> f isa SubArray, factors)
    parents = unique(parent.(factors))
    @test length(parents) == 1
    @test only(parents) isa MB.TrigEvalMatrix
end

function test_bridge_factor_matches_transformation_row_monomial()
    inner, set = _bridge_into_inner(MB.Monomial)
    declared = _added_constraint_set_type(
        MOI.Bridges.Variable.concrete_bridge_type(
            SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
            typeof(set),
        ),
    )
    ci = first(
        MOI.get(
            inner,
            MOI.ListOfConstraintIndices{MOI.VectorOfVariables,declared}(),
        ),
    )
    s = MOI.get(inner, MOI.ConstraintSet(), ci)
    U = MB.transformation_to(only(set.gram_bases), set.basis)
    for j in eachindex(set.basis)
        # Each factor is the j-th *row* of `U`, accessed lazily through a
        # `SubArray` (no materialisation).
        @test s.vectors[j].matrix.factor == view(U, j, :)
    end
end

# 6 ─────────────────────────────────────────────────────────────────────────
# `PolyJuMP.bridges(::Type{<:LRO.SetDotProducts{...}})` — the dispatch hook
# in `src/variables.jl` was retyped to accept *any* `F<:AbstractVector{T}`
# factor type (instead of hardcoded `Vector{T}`). Cover both the previous
# `Vector{T}` and the new `SubArray` shape; both must hit the same bridge
# list so the bridge graph routes them identically.
# ─────────────────────────────────────────────────────────────────────────

function _polyjump_bridges_for(SetType)
    return PolyJuMP.bridges(SetType)
end

function test_polyjump_bridges_dispatches_for_vector_factor()
    T = Float64
    V = LRO.TriangleVectorization{T,LRO.Factorization{T,Vector{T},Array{T,0}}}
    SetType = LRO.SetDotProducts{
        LRO.WITHOUT_SET,
        MOI.PositiveSemidefiniteConeTriangle,
        V,
        Vector{V},
    }
    bridges = _polyjump_bridges_for(SetType)
    @test bridges isa Vector
    @test !isempty(bridges)
    # The classic chain through `ToPositiveBridge` + `AppendSetBridge`.
    @test any(b -> b[1] === LRO.Bridges.Variable.ToPositiveBridge, bridges)
end

function test_polyjump_bridges_dispatches_for_subarray_factor()
    T = Float64
    U = Matrix{Float64}(undef, 3, 2)
    F = typeof(view(U, 1, :))
    V = LRO.TriangleVectorization{T,LRO.Factorization{T,F,Array{T,0}}}
    SetType = LRO.SetDotProducts{
        LRO.WITHOUT_SET,
        MOI.PositiveSemidefiniteConeTriangle,
        V,
        Vector{V},
    }
    bridges = _polyjump_bridges_for(SetType)
    @test bridges isa Vector
    @test !isempty(bridges)
    @test any(b -> b[1] === LRO.Bridges.Variable.ToPositiveBridge, bridges)
end

function test_polyjump_bridges_dispatches_for_trig_eval_subarray_factor()
    T = Float64
    # Points must stay in `[-1, 1]` because `Trigonometric.recurrence_eval`
    # computes `sqrt(1 - cos²θ)` for sin terms. The `_weighted_sos_cone`
    # helper above keeps points in that range; we mirror its layout here.
    _, lag, gram, _ = _weighted_sos_cone(MB.Trigonometric; gram_degree = 4, n_pts = 9)
    U = MB.transformation_to(gram, lag)
    F = typeof(view(U, 1, :))
    @test F <: SubArray
    @test eltype(parent(view(U, 1, :))) === T
    V = LRO.TriangleVectorization{T,LRO.Factorization{T,F,Array{T,0}}}
    SetType = LRO.SetDotProducts{
        LRO.WITHOUT_SET,
        MOI.PositiveSemidefiniteConeTriangle,
        V,
        Vector{V},
    }
    bridges = _polyjump_bridges_for(SetType)
    @test bridges isa Vector
    @test !isempty(bridges)
end

# ─────────────────────────────────────────────────────────────────────────

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

end

TestVariableLowRankParam.runtests()
