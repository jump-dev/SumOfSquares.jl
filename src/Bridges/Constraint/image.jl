"""
    ImageBridge{T,F,G,MT,MVT,CT} <: Bridges.Constraint.AbstractBridge

`ImageBridge` implements a reformulation from `SOSPolynomialSet{SemialgebraicSets.FullSpace}`
into [`MOI.PositiveSemidefiniteConeTriangle`](@ref).

Let `Σ` be the SOS cone of polynomials of degree 2d and `S` be the PSD cone.
There is a linear relation `Σ = A(S)`.
The linear relation reads: `p` belongs to `Σ` iff there exists `q` in `S` such that `A(q) = p`.
This allows defining a variable bridge that would create variables `p` and substitute `A(q)` for `p` but this is not the purpose of this bridge.
This bridge exploit the following alternative read:
`p` belongs to `Σ` iff there exists `q` in `S` such that ``q in A^{-1}(p)`` where `A^{-1}` is the preimage of `p`.
This preimage can be obtained as `A^\\dagger p + \\mathrm{ker}(A)` where `A^\\dagger` is the pseudo-inverse of `A`.
It turns out that for polynomial bases indexed by monomials, `A` is close to row echelon form so
`A^\\dagger` and `\\mathrm{ker}(A)` can easily be obtained.

This is best described in an example.
Consider the SOS constraint for the polynomial `p = 2x^4 + 2x^3 * y - x^2 * y^2 + 5y^4`
with gram basis `b = [x^2, y^2, x * y]` of [Parrilo2003; Example 6.1](@cite).
The product `b * b'` is
```math
\\begin{bmatrix}
x^4 & x^2 y^2 & x^3 y\\
x^2 y^2 & y^4 & x y^3\\
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
2 & -\\lambda & 1\\
-\\lambda & 5 & 0\\
1 & 0 & 2\\lambda - 1
\\end{bmatrix}
```

## Source node

`ImageBridge` supports:

  * `H` in `SOSPolynomialSet{SemialgebraicSets.FullSpace}`

## Target nodes

`ImageBridge` creates one of the following, depending on the length of the gram basis:

  * `F` in `MOI.PositiveSemidefiniteConeTriangle`, for gram basis of length at least 3
  * `F` in [`PositiveSemidefinite2x2ConeTriangle`](@ref), for gram basis of length 2
  * `F` in `MOI.Nonnegatives`, for gram basis of length 1
  * `F` in `EmptyCone`, for empty gram basis

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
    k = 0
    found = Dict{eltype(set.basis.monomials),Int}()
    first = Union{Nothing,Int}[nothing for _ in eachindex(scalars)]
    variables = MOI.VariableIndex[]
    constraints = MOI.ConstraintIndex{F}[]
    for (gram_basis, weight) in zip(set.gram_bases, set.weights)
        cone = SOS.matrix_cone(M, length(gram_basis))
        f = MOI.Utilities.zero_with_output_dimension(F, MOI.dimension(cone))
        for j in eachindex(gram_basis.monomials)
            for i in 1:j
                k += 1
                mono = gram_basis.monomials[i] * gram_basis.monomials[j]
                is_diag = i == j
                if haskey(found, mono)
                    var = MOI.add_variable(model)
                    push!(variables, var)
                    is_diag_found =
                        MOI.Utilities.is_diagonal_vectorized_index(found[mono])
                    if is_diag == is_diag_found
                        MOI.Utilities.operate_output_index!(
                            +,
                            T,
                            found[mono],
                            f,
                            var,
                        )
                    else
                        coef = is_diag ? inv(T(2)) : T(2)
                        MOI.Utilities.operate_output_index!(
                            +,
                            T,
                            found[mono],
                            f,
                            coef * var,
                        )
                    end
                    MOI.Utilities.operate_output_index!(-, T, k, f, var)
                else
                    found[mono] = k
                    t = MP.searchsortedfirst(set.basis.monomials, mono)
                    if t in eachindex(set.basis.monomials) && set.basis.monomials[t] == mono
                        first[t] = k
                        if is_diag
                            MOI.Utilities.operate_output_index!(
                                +,
                                T,
                                k,
                                f,
                                scalars[t],
                            )
                        else
                            MOI.Utilities.operate_output_index!(
                                +,
                                T,
                                k,
                                f,
                                inv(T(2)) * scalars[t],
                            )
                        end
                    end
                end
            end
        end
        push!(constraints, MOI.add_constraint(model, f, cone))
    end
    if any(isnothing, first)
        z = findall(isnothing, first)
        zero_constraint = MOI.add_constraint(
            model,
            MOI.Utilities.vectorize(scalars[z]),
            MOI.Zeros(length(z)),
        )
    else
        zero_constraint = nothing
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
    ::Type{<:SOS.WeightedSOSCone{M,<:MB.MonomialBasis,<:MB.MonomialBasis}},
) where {T,M}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:ImageBridge})
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:ImageBridge{T,F,G}},
) where {T,F,G}
    return Tuple{Type,Type}[
        (F, MOI.PositiveSemidefiniteConeTriangle),
        (G, MOI.Zeros),
    ] # TODO
end

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
        ci isa MOI.ConstraintIndex{F,S}
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
        ci for ci in bridge.constraints
        if ci isa MOI.ConstraintIndex{F,S}
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

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintFunction, bridge::ImageBridge{T}) where {T}
    if !isnothing(bridge.zero_constraint)
        z = MOI.Utilities.eachscalar(MOI.get(model, attr, bridge.zero_constraint))
    end
    funcs = MOI.Utilities.eachscalar.(MOI.get.(model, MOI.ConstraintFunction(), bridge.constraints))
    z_idx = 0
    return MOI.Utilities.vectorize(map(eachindex(bridge.first)) do i
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
    end)
end

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintPrimal, ::ImageBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintDual,PolyJuMP.MomentsAttribute},
    bridge::ImageBridge{T},
) where {T}
    dual =
        MOI.get(model, MOI.ConstraintDual(attr.result_index), bridge.constraint)
    output = similar(dual, length(bridge.set.monomials))
    for i in eachindex(bridge.set.monomials)
        output[i] = dual[bridge.first[i]]
    end
    return output
end

function MOI.get(
    ::MOI.ModelLike,
    ::SOS.CertificateBasis,
    bridge::ImageBridge,
)
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
    return SOS.build_gram_matrix(
        q,
        bridge.gram_basis,
        M,
        T,
    )
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
