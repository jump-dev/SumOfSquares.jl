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

### Non-unit weight

When the weight `w` of the gram basis is not one, each entry of `b * b'`
contributes to a polynomial coefficient with an additional factor `w`.
Concretely, every `(i, j)` entry of the gram matrix gets divided by `w`
(and by the usual extra `2` for off-diagonal entries), while the slack
adjustment ratio `(factor / factor_anchor)` is unchanged since `w` cancels.

Take the same polynomial `p` and gram basis as above, but now write the
SOS constraint with weight `w = 1 + x` and `q = q_0 + q_1 x`, where each `q_i`
is a coefficient of the gram matrix.
For brevity, denote the polynomial coefficients
`p = a_0 + a_1 x + a_2 x^2 + a_3 x^3` with gram basis `b = [1, x]`.
The product `b * b'` is
```math
\\begin{bmatrix}
1 & x\\\\
x & x^2
\\end{bmatrix}
```
and `(1 + x) * b * b'` distributes one extra power of `x`:
the `(1, 1)` entry contributes both to the `1` and the `x` coefficient of `p`,
the `(1, 2)` (and `(2, 1)`) entry contributes to `x` and `x^2`, and the
`(2, 2)` entry contributes to `x^2` and `x^3`.
Anchoring each monomial of `p` to the first gram entry that produces it
(while accumulating the contributions of every earlier entry) gives:
- `q_{1,1} = a_0`;
- `q_{1,2}` is anchored on the `x` coefficient and absorbs `q_{1,1}`'s
  contribution, yielding `q_{1,2} = (a_1 - a_0) / 2`;
- `q_{2,2}` is anchored on the `x^2` coefficient and absorbs the contribution
  of `q_{1,2}`, yielding `q_{2,2} = a_2 - a_1 + a_0`.
There is no entry left to anchor the `x^3` coefficient, so its equation
`a_3 = q_{2,2}` becomes the zero constraint
`a_0 - a_1 + a_2 - a_3 = 0`.
The gram matrix is therefore:
```math
\\begin{bmatrix}
a_0 & \\frac{a_1 - a_0}{2}\\\\
\\frac{a_1 - a_0}{2} & a_0 - a_1 + a_2
\\end{bmatrix}
```
with the additional zero constraint `a_0 - a_1 + a_2 - a_3 = 0`.

### Multiple weights and gram bases

The same logic extends to a `WeightedSOSCone` with several gram bases
`b_i` and weights `w_i`: the bridge walks through every `(i, j, k)` entry
of every gram matrix and treats each weight monomial as an additional
contribution to the same anchoring/slack bookkeeping. Anchors and slack
variables can therefore span across different gram matrices.

Consider for instance Putinar-style decomposition
`p = σ_0 + x * σ_1` with
`p = a_0 + a_1 x + a_2 x^2`,
weights `w_0 = 1`, `w_1 = x`,
and gram bases `b_0 = [1, x]` (a 2 × 2 gram matrix `Q`)
and `b_1 = [1]` (a 1 × 1 gram matrix `R`).

Walking through `b_0`'s entries first (weight `1`):
- `Q_{1,1}` anchors `1`, so `Q_{1,1} = a_0`;
- `Q_{1,2}` anchors `x` and `Q_{1,2} = a_1 / 2`;
- `Q_{2,2}` anchors `x^2` and `Q_{2,2} = a_2`.

Then `R_{1,1}` with weight `x` produces the monomial
`x · 1 · 1 = x`, which is already anchored by `Q_{1,2}` (off-diagonal,
factor `2`). Introducing a slack variable `λ` to balance, with adjustment
ratio `1 / 2`:
`R_{1,1} = -λ` and `Q_{1,2}` is updated to `(a_1 + λ) / 2`.
The gram matrices are therefore
```math
Q = \\begin{bmatrix}
a_0 & \\frac{a_1 + \\lambda}{2}\\\\
\\frac{a_1 + \\lambda}{2} & a_2
\\end{bmatrix},
\\qquad
R = \\begin{bmatrix}-\\lambda\\end{bmatrix}.
```
The slack `λ` is constrained only by `Q ⪰ 0` and `R ⪰ 0`,
i.e. by the PSD constraints on the gram matrices.

### Sub-optimal greedy: slacks that could be avoided

The order in which the contributions of one entry are processed matters.
The bridge walks each entry's weight monomials in their stored order and
commits to either anchoring the entry to the first unanchored target it
sees, or to a slack variable if that first target is already anchored.
That decision can introduce a slack variable that a different order would
have avoided.

Consider weight `w = 1 + x` with gram basis `[1, x, x^2]` and basis
`{1, x, x^2, x^3, x^4, x^5}`. The ``6`` PSD entries are
```text
k = 1: (1, 1)     contributes to 1   and x
k = 2: (1, x)     contributes to x   and x^2
k = 3: (x, x)     contributes to x^2 and x^3
k = 4: (1, x^2)   contributes to x^2 and x^3
k = 5: (x, x^2)   contributes to x^3 and x^4
k = 6: (x^2, x^2) contributes to x^4 and x^5
```
The greedy walks weight monomials in the order `[1, x]`, so for entry
``k = 4`` it sees the contribution to `x^2` first. By then `x^2` is
already anchored at ``k = 3``, so the bridge introduces a slack `λ`,
records the second contribution (`x^3`, currently unanchored) into
`pending`, and the trail of consequences cascades into 5 anchors and 1
slack — with `x^5` left to a zero constraint.

A bipartite matching between gram entries (left) and basis monomials
(right) — an edge `(k, t)` for each weight-monomial contribution from
``k`` to a target ``t`` in the basis — admits a perfect matching of size
6 here:
``1 \\to 1,\\ 2 \\to x,\\ 3 \\to x^2,\\ 4 \\to x^3,\\ 5 \\to x^4,\\ 6 \\to x^5``.
Anchoring along that matching would use 6 anchors and 0 slack. Even just
swapping the two contributions of entry ``k = 4`` (and similarly for
``k = 5`` and ``k = 6``) — i.e. preferring an unanchored target whenever
one is available, and falling back to subtracting from the existing
anchor for the other contribution — recovers the same outcome. The
underlying problem is unweighted bipartite max-cardinality matching, for
which Hopcroft–Karp or simple augmenting paths gives the optimum; the
weighted Hungarian assignment is not needed.

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
    # `first[t] = (gram_basis_idx, k_local)` identifies which entry in which
    # gram basis anchors the `t`-th basis monomial.
    first::Vector{Union{Nothing,Tuple{Int,Int}}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ImageBridge{T,F,G,M}},
    model::MOI.ModelLike,
    g::MOI.AbstractVectorFunction,
    set::SOS.WeightedSOSCone{M},
) where {T,F,G,M}
    @assert MOI.output_dimension(g) == length(set.basis)
    scalars = MOI.Utilities.scalarize(g)
    # `found[mono] = (i_basis, k, factor)` records that the monomial `mono`
    # is anchored by the `k`-th entry of the `i_basis`-th gram matrix, with
    # the recorded factor.
    found = Dict{eltype(set.basis),Tuple{Int,Int,T}}()
    # `pending[t]` lists `(i_basis, k, factor)` contributions to the `t`-th
    # basis monomial that is not yet anchored; the contributions will be
    # absorbed when the monomial gets anchored, or otherwise turned into a
    # zero constraint.
    pending = Dict{Int,Vector{Tuple{Int,Int,T}}}()
    first = Union{Nothing,Tuple{Int,Int}}[nothing for _ in eachindex(scalars)]
    variables = MOI.VariableIndex[]
    # The gram matrices are built into `fs` first, then handed to MOI at
    # the end: when an entry from a later gram basis slacks against the
    # anchor of an earlier basis, the earlier `f` still has to be mutable.
    fs = F[]
    for gram_basis in set.gram_bases
        cone = SOS.matrix_cone(M, length(gram_basis))
        push!(
            fs,
            MOI.Utilities.zero_with_output_dimension(F, MOI.dimension(cone)),
        )
    end
    for (i_basis, (gram_basis, weight)) in
        enumerate(zip(set.gram_bases, set.weights))
        weight_basis = SA.basis(weight)
        weight_coeffs = SA.coeffs(weight)
        f = fs[i_basis]
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
                        i_a, k_a, f_a = found[mono]
                        f_a_func = fs[i_a]
                        if entry_anchored
                            f_k = MOI.Utilities.eachscalar(f)[k]
                            MOI.Utilities.operate_output_index!(
                                -,
                                T,
                                k_a,
                                f_a_func,
                                (f_kw / f_a) * f_k,
                            )
                        else
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
                                f_a_func,
                                (f_kw / f_a) * slack_var,
                            )
                        end
                    elseif mono in set.basis
                        t = set.basis[mono]
                        if !entry_anchored && slack_var === nothing
                            found[mono] = (i_basis, k, f_kw)
                            first[t] = (i_basis, k)
                            MOI.Utilities.operate_output_index!(
                                +,
                                T,
                                k,
                                f,
                                inv(f_kw) * scalars[t],
                            )
                            if haskey(pending, t)
                                for (i_p, k_p, f_p) in pending[t]
                                    f_kp =
                                        MOI.Utilities.eachscalar(fs[i_p])[k_p]
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
                            push!(
                                get!(pending, t, Tuple{Int,Int,T}[]),
                                (i_basis, k, f_kw),
                            )
                        end
                    end
                end
            end
        end
    end
    constraints = MOI.ConstraintIndex{F}[]
    for (i_basis, gram_basis) in enumerate(set.gram_bases)
        cone = SOS.matrix_cone(M, length(gram_basis))
        push!(constraints, MOI.add_constraint(model, fs[i_basis], cone))
    end
    # Build zero constraints for unanchored basis monomials.
    z = findall(isnothing, first)
    zero_terms = if isempty(z)
        empty(scalars)
    else
        map(z) do t
            term = scalars[t]
            if haskey(pending, t)
                for (i_p, k_p, f_p) in pending[t]
                    f_kp = MOI.Utilities.eachscalar(fs[i_p])[k_p]
                    term = MA.operate!!(-, term, f_p * f_kp)
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
    ::Type{<:SOS.WeightedSOSCone{M,B}},
) where {T,M,B}
    # `LagrangeBasis` dispatches to `Variable.LowRankBridge`; the image-form
    # bridge cannot iterate over its samples (`bridge_constraint` walks the
    # source basis monomial-by-monomial), so excluding it here both fixes the
    # runtime error and lets the cost graph correctly favour `LowRankBridge`.
    return !(B <: MB.LagrangeBasis)
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
                i_basis, k = bridge.first[i]
                f = MOI.Utilities.filter_variables(
                    !Base.Fix2(in, bridge.variables),
                    funcs[i_basis][k],
                )
                if !MOI.Utilities.is_diagonal_vectorized_index(k)
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
