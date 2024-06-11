module Certificate

import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariatePolynomials as MP
import MultivariateBases as MB
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

# For get
function cone end
function gram_basis end
# FIXME currently, this returns `MB.MonomialBasis` instead of `MB.MonomialBasis{MT, MVT}`
function gram_basis_type end
function reduced_polynomial end
function ideal_certificate end
function preprocessed_domain end

# Only for PreorderIndex
function preorder_indices end
function generator end
function multiplier_basis end
function multiplier_basis_type end

abstract type AbstractCertificate end

function within_total_bounds(mono::MP.AbstractMonomial, bounds::DegreeBounds)
    return bounds.mindegree <= MP.degree(mono) <= bounds.maxdegree
end

function within_variablewise_bounds(
    mono::MP.AbstractMonomial,
    bounds::DegreeBounds,
)
    return MP.divides(bounds.variablewise_mindegree, mono) &&
           MP.divides(mono, bounds.variablewise_maxdegree)
end

function within_bounds(mono, bounds)
    return within_total_bounds(mono, bounds) &&
           within_variablewise_bounds(mono, bounds)
end

function maxdegree_gram_basis(
    ::MB.FullBasis{B},
    bounds::DegreeBounds,
) where {B<:MB.AbstractMonomial}
    variables = MP.variables(bounds.variablewise_maxdegree)
    return MB.SubBasis{B}(
        MP.monomials(
            variables,
            bounds.mindegree:bounds.maxdegree,
            Base.Fix2(within_variablewise_bounds, bounds),
        ),
    )
end

function maxdegree_gram_basis(basis::SA.AbstractBasis, ::Nothing)
    return MB.empty_basis(MB.explicit_basis_type(typeof(basis)))
end

function maxdegree_gram_basis(basis::SA.AbstractBasis, bounds::DegreeBounds)
    # TODO use bounds here too
    variables = MP.variables(bounds.variablewise_maxdegree)
    return MB.maxdegree_basis(basis, variables, bounds.maxdegree)
end
function maxdegree_gram_basis(
    basis::SA.AbstractBasis,
    variables,
    maxdegree::Int,
)
    return MB.maxdegree_basis(basis, variables, fld(maxdegree, 2))
end

include("ideal.jl")
include("preorder.jl")

include("Sparsity/Sparsity.jl")
using .Sparsity: SignSymmetry, ChordalCompletion, ClusterCompletion
export Sparsity, SignSymmetry, ChordalCompletion, ClusterCompletion

include("Symmetry/Symmetry.jl")
export Symmetry

end
