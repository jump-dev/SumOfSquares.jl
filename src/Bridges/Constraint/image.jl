"""
    ImageBridge{T,F,G,MT,MVT,CT} <: Bridges.Constraint.AbstractBridge

`ImageBridge` implements a reformulation from `SOSPolynomialSet{SemialgebraicSets.FullSpace}`
into the positive semidefinite cone.

Let `Σ` be the SOS cone of polynomials of degree 2d and `S` be the PSD cone.
There is a linear relation `Σ = A(S)`.
The linear relation reads: `p` belongs to `Σ` iff there exists `q` in `S` such that `A(q) = p`.
This allows defining a variable bridge that would create variables `p` and substitute `A(q)` for `p` but this is not the purpose of this bridge.
This bridge exploit the following alternative read:
`p` belongs to `Σ` iff there exists `q` in `S` such that ``q \\in A^{-1}(p)`` where ``A^{-1}`` is the preimage of `p`.
This preimage can be obtained as ``A^\\dagger p + \\mathrm{ker}(A)`` where ``A^\\dagger`` is the pseudo-inverse of `A`.
It turns out that for polynomial bases indexed by monomials, `A` is close to row echelon form so
``A^\\dagger`` and ``\\mathrm{ker}(A)`` can easily be obtained.

This is best described in an example.
Consider the SOS constraint for the polynomial `p = 2x^4 + 2x^3 * y - x^2 * y^2 + 5y^4`
with gram basis `b = [x^2, y^2, x * y]` of [Parrilo2003; Example 6.1](@cite).
The product `b * b'` is
```math
\\begin{bmatrix}
x^4 & x^2 y^2 & x^3 y\\\\
x^2 y^2 & y^4 & x y^3\\\\
x^3 y & x y^3 & x^2 y^2
\\end{bmatrix}
```
Except for the entries `(1, 2)` and `(3, 3)`, all entries are unique so the value of the
corresponding gram matrix is given by the corresponding coefficient of `p`.
For the entries `(1, 2)` and `(3, 3)`, we let `-λ` be the `(1, 2)` and `(3, 3)` entries
and we add `2λ` for the `(3, 3)` entry so as to parametrize entries for which the sum is
the corresponding coefficient in `p`, i.e., `-1`.
The gram matrix is therefore:
```math
\\begin{bmatrix}
2 & -\\lambda & 1\\\\
-\\lambda & 5 & 0\\\\
1 & 0 & 2\\lambda - 1
\\end{bmatrix}
```

## Source node

`ImageBridge` supports:

  * `H` in `SOSPolynomialSet{SemialgebraicSets.FullSpace}`

## Target nodes

`ImageBridge` creates one of the following, depending on the length of the gram basis:

  * `F` in `MOI.PositiveSemidefiniteConeTriangle`, for gram basis of length at least 3
  * `F` in [`SumOfSquares.PositiveSemidefinite2x2ConeTriangle`](@ref), for gram basis of length 2
  * `F` in `MOI.Nonnegatives`, for gram basis of length 1
  * `F` in [`SumOfSquares.EmptyCone`](@ref), for empty gram basis

in addition to

  * a constraint `G` in `MOI.Zeros` in case there is a monomial in `s.monomials`
    that cannot be obtained as product of elements in a gram basis.
"""
struct ImageBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    G<:MOI.AbstractVectorFunction,
    M,
} <: MOI.Bridges.Constraint.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    constraints::Vector{MOI.ConstraintIndex{F}}
    zero_constraint::Union{Nothing,MOI.ConstraintIndex{G,MOI.Zeros}}
    set::SOS.WeightedSOSCone{M}
    first::Vector{Union{Nothing,Int}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ImageBridge{T,F,G,M}},
    model::MOI.ModelLike,
    g::MOI.AbstractVectorFunction,
    set::SOS.WeightedSOSCone{M},
) where {T,F,G,M}
    @assert MOI.output_dimension(g) == length(set.basis)
    scalars = MOI.Utilities.scalarize(g)
    # `found[mono] = k` records the index in the current `f` where the
    # monomial `mono` has been anchored, with factor `factor[mono]`.
    found = Dict{eltype(set.basis),Int}()
    factor = Dict{eltype(set.basis),T}()
    # `pending[t]` lists `(k, factor)` contributions to the `t`-th basis
    # monomial that is not yet anchored; the contributions will be absorbed
    # when the monomial gets anchored, or otherwise turned into a zero
    # constraint.
    pending = Dict{Int,Vector{Tuple{Int,T}}}()
    first = Union{Nothing,Int}[nothing for _ in eachindex(scalars)]
    variables = MOI.VariableIndex[]
    constraints = MOI.ConstraintIndex{F}[]
    # `fs[i]` keeps the function passed to `MOI.add_constraint` for the
    # `i`-th gram basis. It is used below to build the residual zero
    # constraint for unanchored basis monomials.
    fs = F[]
    for (gram_basis, weight) in zip(set.gram_bases, set.weights)
        cone = SOS.matrix_cone(M, length(gram_basis))
        f = MOI.Utilities.zero_with_output_dimension(F, MOI.dimension(cone))
        weight_basis = SA.basis(weight)
        weight_coeffs = SA.coeffs(weight)
        k = 0
        for j in eachindex(gram_basis)
            for i in 1:j
                k += 1
                is_diag = i == j
                diag_factor = is_diag ? one(T) : T(2)
                entry_anchored = false
                slack_var = nothing
                for (w_key, w_coef) in
                    zip(SA.keys(weight_coeffs), SA.values(weight_coeffs))
                    mono_w = weight_basis[w_key]
                    mono = mono_w * SA.star(gram_basis[i]) * gram_basis[j]
                    f_kw = w_coef * diag_factor
                    if haskey(found, mono)
                        k_a = found[mono]
                        f_a = factor[mono]
                        if entry_anchored
                            # Entry `k` is anchored by a previous weight term:
                            # subtract its current contribution from the
                            # anchor of `mono`.
                            f_k = MOI.Utilities.eachscalar(f)[k]
                            MOI.Utilities.operate_output_index!(
                                -,
                                T,
                                k_a,
                                f,
                                (f_kw / f_a) * f_k,
                            )
                        else
                            # Entry `k` is free so far: introduce a slack
                            # variable representing its value (or extend the
                            # one already introduced by an earlier weight
                            # term) and adjust the anchor.
                            if slack_var === nothing
                                slack_var = MOI.add_variable(model)
                                push!(variables, slack_var)
                                MOI.Utilities.operate_output_index!(
                                    -,
                                    T,
                                    k,
                                    f,
                                    slack_var,
                                )
                            end
                            MOI.Utilities.operate_output_index!(
                                +,
                                T,
                                k_a,
                                f,
                                (f_kw / f_a) * slack_var,
                            )
                        end
                    elseif mono in set.basis
                        t = set.basis[mono]
                        if !entry_anchored && slack_var === nothing
                            # Anchor `mono` at entry `k` with factor `f_kw`.
                            found[mono] = k
                            factor[mono] = f_kw
                            first[t] = k
                            MOI.Utilities.operate_output_index!(
                                +,
                                T,
                                k,
                                f,
                                inv(f_kw) * scalars[t],
                            )
                            if haskey(pending, t)
                                for (k_p, f_p) in pending[t]
                                    f_kp = MOI.Utilities.eachscalar(f)[k_p]
                                    MOI.Utilities.operate_output_index!(
                                        -,
                                        T,
                                        k,
                                        f,
                                        (f_p / f_kw) * f_kp,
                                    )
                                end
                                delete!(pending, t)
                            end
                            entry_anchored = true
                        else
                            # Entry `k` is busy (already anchored or slack):
                            # record the contribution for when `mono` gets
                            # anchored or for a zero constraint.
                            push!(
                                get!(pending, t, Tuple{Int,T}[]),
                                (k, f_kw),
                            )
                        end
                    end
                end
            end
        end
        push!(fs, f)
        push!(constraints, MOI.add_constraint(model, f, cone))
    end
    # Build zero constraints for unanchored basis monomials. For each such
    # monomial, the polynomial coefficient must equal the sum of the pending
    # contributions; if no contributions exist at all, the coefficient must
    # simply be zero.
    z = findall(isnothing, first)
    zero_terms = if isempty(z)
        empty(scalars)
    else
        # The per-gram-basis `f`s above were built with a local `k`. With a
        # single gram basis (the typical case) the global pending index is
        # the local one. Supporting multiple gram bases would require
        # tracking the originating basis for each pending contribution.
        @assert length(fs) == 1
        f_local = MOI.Utilities.eachscalar(fs[1])
        map(z) do t
            term = scalars[t]
            if haskey(pending, t)
                for (k_p, f_p) in pending[t]
                    term = MA.operate!!(-, term, f_p * f_local[k_p])
                end
            end
            return term
        end
    end
    zero_constraint = if isempty(zero_terms)
        nothing
    else
        MOI.add_constraint(
            model,
            MOI.Utilities.vectorize(zero_terms),
            MOI.Zeros(length(zero_terms)),
        )
    end
    return ImageBridge{T,F,G,M}(
        variables,
        constraints,
        zero_constraint,
        set,
        first,
    )
end

function MOI.supports_constraint(
    ::Type{ImageBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.WeightedSOSCone},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:ImageBridge})
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:ImageBridge{T,F,G,M}},
) where {T,F,G,M}
    return Tuple{Type,Type}[
        (F, typeof(SOS.matrix_cone(M, 0))),
        (F, typeof(SOS.matrix_cone(M, 1))),
        (F, typeof(SOS.matrix_cone(M, 2))),
        (F, typeof(SOS.matrix_cone(M, 3))),
        (G, MOI.Zeros),
    ]
end

# The default cost of `1` makes `ImageBridge` win against `Variable.KernelBridge`
# for solvers that natively support PSD as constrained variables (e.g. Mosek's
# `barvar` matrix variables), where `KernelBridge` would be more natural.
# Bumping the cost by `1` shifts the balance: for solvers that only support PSD
# as a constraint (e.g. Clarabel), the `KernelBridge` fallback path is still
# strictly more expensive, so `ImageBridge` remains the cheapest path; for
# solvers that support PSD as constrained variables, `KernelBridge` now wins.
MOI.Bridges.bridging_cost(::Type{<:ImageBridge}) = 2.0

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ImageBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.WeightedSOSCone{M}},
) where {T,M}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    F = MOI.Utilities.promote_operation(-, T, G, MOI.VectorOfVariables)
    return ImageBridge{T,F,G,M}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::ImageBridge, ::MOI.NumberOfVariables)
    return length(bridge.variables)
end
function MOI.get(bridge::ImageBridge, ::MOI.ListOfVariableIndices)
    return bridge.variables
end
function MOI.get(
    bridge::ImageBridge{T,F,G},
    ::MOI.NumberOfConstraints{F,S},
) where {T,F,G,S}
    return count(bridge.constraints) do ci
        return ci isa MOI.ConstraintIndex{F,S}
    end
end
function MOI.get(
    bridge::ImageBridge{T,F,G},
    ::MOI.NumberOfConstraints{G,MOI.Zeros},
) where {T,F,G}
    return isnothing(bridge.zero_constraint) ? 0 : 1
end

function MOI.get(
    bridge::ImageBridge{T,F,G},
    ::MOI.ListOfConstraintIndices{F,S},
) where {T,F,G,S}
    return MOI.ConstraintIndex{F,S}[
        ci for ci in bridge.constraints if ci isa MOI.ConstraintIndex{F,S}
    ]
end

function MOI.get(
    bridge::ImageBridge{T,F,G},
    ::MOI.ListOfConstraintIndices{G,MOI.Zeros},
) where {T,F,G}
    if isnothing(bridge.zero_constraint)
        return MOI.ConstraintIndex{G,MOI.Zeros}[]
    else
        return [bridge.zero_constraint]
    end
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::ImageBridge)
    if !isnothing(bridge.zero_constraint)
        MOI.delete(model, bridge.zero_constraint)
    end
    MOI.delete(model, bridge.constraints)
    if !isempty(bridge.variables)
        MOI.delete(model, bridge.variables)
    end
    return
end

# Attributes, Bridge acting as a constraint
function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, bridge::ImageBridge)
    return bridge.set
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintFunction,
    bridge::ImageBridge{T},
) where {T}
    # Recovering the original function from the bridged constraints requires
    # inverting the affine map applied per gram entry. For weights that are
    # not the constant one, this map is no longer simply identity-with-2
    # scaling, so we declare unbridging unsupported instead of returning a
    # wrong function.
    if !all(isone, bridge.set.weights)
        throw(
            MOI.GetAttributeNotAllowed(
                attr,
                "The `ImageBridge` does not support recovering the original" *
                " function when the `WeightedSOSCone` has non-unit weights.",
            ),
        )
    end
    if !isnothing(bridge.zero_constraint)
        z = MOI.Utilities.eachscalar(
            MOI.get(model, attr, bridge.zero_constraint),
        )
    end
    funcs = MOI.Utilities.eachscalar.(
        MOI.get.(model, MOI.ConstraintFunction(), bridge.constraints),
    )
    z_idx = 0
    return MOI.Utilities.vectorize(
        map(eachindex(bridge.first)) do i
            if isnothing(bridge.first[i])
                z_idx += 1
                return z[z_idx]
            else
                f = MOI.Utilities.filter_variables(
                    !Base.Fix2(in, bridge.variables),
                    funcs[1][bridge.first[i]], # FIXME
                )
                if !MOI.Utilities.is_diagonal_vectorized_index(bridge.first[i])
                    f = T(2) * f
                end
                return f
            end
        end,
    )
end

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintPrimal, ::ImageBridge)
    return throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintDual,PolyJuMP.MomentsAttribute},
    bridge::ImageBridge{T},
) where {T}
    dual =
        MOI.get(model, MOI.ConstraintDual(attr.result_index), bridge.constraint)
    output = similar(dual, length(bridge.set.basis))
    for i in eachindex(bridge.set.basis)
        output[i] = dual[bridge.first[i]]
    end
    return output
end

function MOI.get(::MOI.ModelLike, ::SOS.CertificateBasis, bridge::ImageBridge)
    return bridge.gram_basis
end

function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.GramMatrixAttribute,
    bridge::ImageBridge{T,F,G,M},
) where {T,F,G,M}
    # TODO there are several ones
    q = MOI.get(
        model,
        MOI.ConstraintPrimal(attr.result_index),
        bridge.constraint,
    )
    return SOS.build_gram_matrix(q, bridge.gram_basis, M, T)
end
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.MomentMatrixAttribute,
    bridge::ImageBridge,
)
    # TODO there are several ones
    return SOS.build_moment_matrix(
        MOI.get(
            model,
            MOI.ConstraintDual(attr.result_index),
            bridge.constraint,
        ),
        bridge.gram_basis,
    )
end
