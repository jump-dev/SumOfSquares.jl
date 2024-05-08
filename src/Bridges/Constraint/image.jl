"""
    GeometricBridge{T,F,MT,MVT,CT} <: Bridges.Constraint.AbstractBridge

`GeometricBridge` implements a reformulation from `SOSPolynomialSet{SemialgebraicSets.FullSpace}`
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

`GeometricBridge` supports:

  * `H` in `SOSPolynomialSet{SemialgebraicSets.FullSpace}`

## Target nodes

`GeometricBridge` creates one of the following, depending on the length of the gram basis:

  * `F` in `MOI.PositiveSemidefiniteConeTriangle`, for gram basis of length at least 3
  * `F` in [`PositiveSemidefinite2x2ConeTriangle`](@ref), for gram basis of length 2
  * `F` in `MOI.Nonnegatives`, for gram basis of length 1
  * `F` in `EmptyCone`, for empty gram basis
"""
struct GeometricBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
    CT<:SOS.Certificate.AbstractIdealCertificate,
} <: MOI.Bridges.Constraint.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    constraint::MOI.ConstraintIndex{F}
    set::SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}
    gram_basis::MB.MonomialBasis{MT,MVT}
    first::Vector{Int}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{GeometricBridge{T,F,MT,MVT,CT}},
    model::MOI.ModelLike,
    g::MOI.AbstractVectorFunction,
    s::SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace},
) where {T,F,CT,MT,MVT}
    @assert MOI.output_dimension(g) == length(s.monomials)
    scalars = MOI.Utilities.scalarize(g)
    p = MP.polynomial(scalars, copy(s.monomials))
    gram_basis = SOS.Certificate.gram_basis(
        s.certificate,
        SOS.Certificate.with_variables(p, s.domain),
    )
    MCT = SOS.matrix_cone_type(CT)
    set = SOS.matrix_cone(MCT, length(gram_basis))
    f = MOI.Utilities.zero_with_output_dimension(F, MOI.dimension(set))
    k = 0
    found = Dict{eltype(gram_basis.monomials),Int}()
    first = Union{Nothing,Int}[nothing for _ in eachindex(scalars)]
    variables = MOI.VariableIndex[]
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
                t = MP.searchsortedfirst(s.monomials, mono)
                if t in eachindex(s.monomials) && s.monomials[t] == mono
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
    for t in eachindex(scalars)
        if isnothing(first[t])
            error("Infeasible")
        end
    end
    constraint = MOI.add_constraint(model, f, set)
    return GeometricBridge{T,F,MT,MVT,CT}(
        variables,
        constraint,
        s,
        gram_basis,
        first,
    )
end

function MOI.supports_constraint(
    ::Type{GeometricBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}},
) where {T,CT,MT,MVT}
    return Certificate.gram_basis_type(CT) <: MB.MonomialBasis
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:GeometricBridge})
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:GeometricBridge{T,F}},
) where {T,F}
    return Tuple{Type,Type}[(F, MOI.PositiveSemidefiniteConeTriangle)] # TODO
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:GeometricBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}},
) where {T,MT,MVT,CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    F = MOI.Utilities.promote_operation(-, T, G, MOI.VectorOfVariables)
    return GeometricBridge{T,F,MT,MVT,CT}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::GeometricBridge, ::MOI.NumberOfVariables)
    return length(bridge.variables)
end
function MOI.get(bridge::GeometricBridge, ::MOI.ListOfVariableIndices)
    return bridge.variables
end
function MOI.get(
    bridge::GeometricBridge{T,F},
    ::MOI.NumberOfConstraints{F,S},
) where {T,F,S}
    return bridge.constraint isa MOI.ConstraintIndex{F,S} ? 1 : 0
end

function MOI.get(
    bridge::GeometricBridge{T,F},
    ::MOI.ListOfConstraintIndices{F,S},
) where {T,F,S}
    if bridge.constraint isa MOI.ConstraintIndex{F,S}
        return [bridge.constraint]
    else
        return MOI.ConstraintIndex{F,S}[]
    end
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::GeometricBridge)
    MOI.delete(model, bridge.constraint)
    if !isempty(bridge.variables)
        MOI.delete(model, bridge.variables)
    end
    return
end

# Attributes, Bridge acting as a constraint
function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, bridge::GeometricBridge)
    return bridge.set
end
function MOI.get(::MOI.ModelLike, ::MOI.ConstraintPrimal, ::GeometricBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintDual,PolyJuMP.MomentsAttribute},
    bridge::GeometricBridge{T},
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
    bridge::GeometricBridge,
)
    return bridge.gram_basis
end

function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.GramMatrixAttribute,
    bridge::GeometricBridge{T,F,MT,MVT,CT},
) where {T,F,MT,MVT,CT}
    q = MOI.get(
        model,
        MOI.ConstraintPrimal(attr.result_index),
        bridge.constraint,
    )
    return SOS.build_gram_matrix(
        q,
        bridge.gram_basis,
        SOS.matrix_cone_type(CT),
        T,
    )
end
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.MomentMatrixAttribute,
    bridge::GeometricBridge,
)
    return SOS.build_moment_matrix(
        MOI.get(
            model,
            MOI.ConstraintDual(attr.result_index),
            bridge.constraint,
        ),
        bridge.gram_basis,
    )
end
