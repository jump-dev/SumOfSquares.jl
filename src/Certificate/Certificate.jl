module Certificate

import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM
using SemialgebraicSets

using SumOfSquares

include("newton_polytope.jl")

# Certificate              denominator * p = ...                                      domain
# Polya       (1 + x1 + ... + xn)^(2d) * p = ∑ λ_α x.^α                  λ_α ∈ R_+    R_+^n
# Putinar                                p = σ_0 + ∑ σ_i g_i             σ_i SOS      g_i ≥ 0
# Krivine                                p = ∑ λ_{αβ} g.^α (1 .- g).^β   λ_{αβ} ≥ 0   0 ≤ g_i ≤ 1
# Schmüdgen                              p = ∑ σ_α g^α                   σ_α SOS      g_i ≥ 0
# Hilbert                            d * p = σ                           d, σ SOS     R^n

struct PreorderIndex
    value::Int
end
struct IdealIndex
    value::Int
end

const Index = Union{PreorderIndex,IdealIndex}

abstract type Attribute end

function cone end
function zero_basis end
function gram_basis end
"""
    gram_weights(certificate, gram_basis, poly, ::Type{T})

Return the (parallel) weights for each gram basis returned by `gram_basis`.
The default is one constant `1` weight (single basis case). Certificates that
return a `Vector` of bases (e.g. `Symmetry.Ideal`) can override this to attach
a different weight to each basis (e.g. `degree(χ)` for symmetry-adapted blocks).
"""
function gram_weights end
function reduced_polynomial end
function ideal_certificate end
function preprocessed_domain end

# Only for PreorderIndex
function preorder_indices end
function generator end
function multiplier_basis end
function multiplier_basis_type end

abstract type AbstractCertificate end

# By default, the gram basis does not depend on the domain
gram_basis(certificate, poly, _) = gram_basis(certificate, poly)

function within_total_bounds(mono::MP.AbstractMonomial, bounds::DegreeBounds)
    return bounds.mindegree <= MP.degree(mono) <= bounds.maxdegree
end

_vec(v::AbstractVector) = v
# For `TypedPolynomials`
_vec(v::Tuple) = MP.variable_union_type(first(v))[v...]

function _divides(a, b)
    # `MP.divides(a, b)` is not implemented yet for noncommutative
    vars = unique!(sort(_vec(MP.variables(a))))
    comm = MP.is_commutative(vars)
    return all(vars) do v
        return _degree(a, v, comm) <= _degree(b, v, comm)
    end
end

function within_variablewise_bounds(
    mono::MP.AbstractMonomial,
    bounds::DegreeBounds,
)
    return _divides(bounds.variablewise_mindegree, mono) &&
           _divides(mono, bounds.variablewise_maxdegree)
end

function within_bounds(mono, bounds)
    return within_total_bounds(mono, bounds) &&
           within_variablewise_bounds(mono, bounds)
end

function maxdegree_gram_basis(
    full::MB.FullBasis{B},
    bounds::DegreeBounds,
) where {B<:MB.AbstractMonomialBasis}
    variables = MP.variables(bounds.variablewise_maxdegree)
    monos = MP.monomials(
        variables,
        (bounds.mindegree):(bounds.maxdegree),
        Base.Fix2(within_variablewise_bounds, bounds),
    )
    sub = MB.SubBasis{B}(monos)
    new_sub, new_full = SA.promote_bases(sub, full)
    @assert new_full === full
    return new_sub
end

function maxdegree_gram_basis(basis::SA.AbstractBasis, ::Nothing)
    return SA.SubBasis(basis, SA.key_type(basis)[])
end

function maxdegree_gram_basis(basis::SA.AbstractBasis, bounds::DegreeBounds)
    # TODO use bounds here too
    @assert MP.variables(basis) == MP.variables(bounds.variablewise_maxdegree)
    return MB.maxdegree_basis(basis, bounds.maxdegree)
end

function maxdegree_gram_basis(
    basis::SA.AbstractBasis,
    variables,
    maxdegree::Int,
)
    if MP.variables(basis) == variables
        return MB.maxdegree_basis(basis, fld(maxdegree, 2))
    else
        return _maxdegree_gram_basis(basis, variables, fld(maxdegree, 2))
    end
end

function _maxdegree_gram_basis(
    full::MB.FullBasis{B},
    variables,
    halfdegree::Int,
) where {B}
    monos = MP.monomials(variables, 0:halfdegree)
    sub = MB.SubBasis{B}(monos)
    new_sub, _ = SA.promote_bases(sub, full)
    return new_sub
end

# Gram basis of `σ` such that the degree of `σ` is at most `maxdegree` for a
# certificate where `σ` only needs to be determined modulo the ideal of
# `domain`. If a Gröbner basis of the ideal is known, each polynomial `q` of
# `σ = ∑ q^2` can be replaced by its remainder `r` modulo the Gröbner basis as
# `q^2 - r^2 = (q + r) * (q - r)` is in the ideal. The degree of `r` is at most
# the degree of `q` so only the standard monomials of degree at most
# `maxdegree / 2` are needed, see [Section 4, PP12].
#
# [PP12] Permenter, Frank and Parrilo, Pablo A.
# *Selecting a monomial basis for sums of squares programming over a quotient
# ring*. 51st IEEE Conference on Decision and Control (2012): 1871-1876.
function maxdegree_gram_basis(basis, variables, maxdegree::Int, _)
    return maxdegree_gram_basis(basis, variables, maxdegree)
end

function maxdegree_gram_basis(
    basis::MB.FullBasis{MB.Monomial},
    variables,
    maxdegree::Int,
    domain::AlgebraicSet,
)
    I = SemialgebraicSets.ideal(domain)
    if iszero(SemialgebraicSets.nequalities(domain)) ||
       !SemialgebraicSets.is_gröbner_basis(I)
        return maxdegree_gram_basis(basis, variables, maxdegree)
    end
    monos =
        SemialgebraicSets.standard_monomials(I, variables, fld(maxdegree, 2))
    sub, _ = SA.promote_bases(MB.SubBasis{MB.Monomial}(monos), basis)
    return sub
end

function maxdegree_gram_basis(
    basis,
    variables,
    maxdegree::Int,
    domain::BasicSemialgebraicSet,
)
    return maxdegree_gram_basis(basis, variables, maxdegree, domain.V)
end

include("ideal.jl")
include("preorder.jl")

include("Sparsity/Sparsity.jl")
using .Sparsity: SignSymmetry, ChordalCompletion, ClusterCompletion
export Sparsity, SignSymmetry, ChordalCompletion, ClusterCompletion

include("Symmetry/Symmetry.jl")
export Symmetry

end
