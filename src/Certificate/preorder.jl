#########################
# Preorder certificates #
#########################

abstract type AbstractPreorderCertificate <: AbstractCertificate end

"""
    struct Putinar{
        MC<:AbstractIdealCertificate,
        IC<:AbstractIdealCertificate,
    } <: AbstractPreorderCertificate
        multipliers_certificate::MC
        ideal_certificate::IC
        maxdegree::Int
    end

The `Putinar` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`g_i(x) >= 0` and `h_i(x) = 0` by exhibiting Sum-of-Squares polynomials `σ_i(x)`
such that `p(x) - ∑ σ_i(x) g_i(x)` is guaranteed to be nonnegativity for all `x`
such that `h_i(x) = 0`.
The polynomials `σ_i(x)` are search over `cone` with a basis of type `basis` such that
the degree of `σ_i(x) g_i(x)` does not exceed `maxdegree`.
"""
struct Putinar{MC<:AbstractIdealCertificate,IC<:AbstractIdealCertificate} <:
       AbstractPreorderCertificate
    multipliers_certificate::MC
    ideal_certificate::IC
    # FIXME `maxdegree` here is needed if `multipliers_certificate` is `Newton`
    #       but duplicate if `multipliers_certificate` is `MaxDegree`
    maxdegree::Int
end

cone(certificate::Putinar) = cone(certificate.multipliers_certificate)

# SemialgebraicSets.FullSpace does not implement MP.variables so we need this
# It would be cleaner to make FullSpace have a list of variables though
struct WithVariables{S,V}
    inner::S
    variables::V
end

function MP.variables(v::WithVariables)
    return v.variables
end
SA.basis(v::WithVariables) = SA.basis(v.inner)
MB.explicit_basis(v::WithVariables) = MB.explicit_basis(v.inner)
_algebra_element(v::WithVariables) = v.inner
_algebra_element(a::SA.AlgebraElement) = a

struct WithFixedBases{S,B}
    inner::S
    bases::Vector{B}
    # `σ_0`'s gram basis as computed by `half_newton_polytope` from `p` and
    # the inequality polynomials `gs` (so it already accounts for the
    # multiplier-times-g terms that show up in the SOS decomposition).
    ideal_basis::B
end

# TODO temporary workaround because SS doesn't support AlgebraElement yet
function with_variables(p::SA.AlgebraElement, ::FullSpace)
    return WithVariables(p, MP.variables(p))
end

function with_variables(p::SA.AlgebraElement, domain)
    _, q = SumOfSquares._promote_bases(domain, p)
    return WithVariables(q, MP.variables(q))
end

# TODO not needed
function with_variables(p::MP.AbstractPolynomialLike, ::FullSpace)
    return WithVariables(p, MP.variables(p))
end

function with_variables(p::MP.AbstractPolynomialLike, domain)
    inner, outer = SumOfSquares._promote_bases(p, domain)
    return WithVariables(inner, MP.variables(outer))
end

function with_variables(domain, p::WithVariables)
    return with_variables(domain, p.inner)
end

function with_variables(domain, p)
    inner, outer = SumOfSquares._promote_bases(domain, p)
    return WithVariables(inner, MP.variables(outer))
end

function with_fixed_basis(
    domain,
    p,
    maxdegree,
    newton::AbstractNewtonPolytopeApproximation,
)
    v = with_variables(domain, p)
    ideal_basis, multiplier_bases = half_newton_polytope(
        _algebra_element(p),
        SemialgebraicSets.inequalities(v.inner),
        v.variables,
        maxdegree,
        newton,
    )
    return WithFixedBases(v.inner, multiplier_bases, ideal_basis)
end

function preprocessed_domain(
    ::Putinar{<:MaxDegree},
    domain::BasicSemialgebraicSet,
    p,
)
    return with_variables(domain, p)
end

function preprocessed_domain(
    certificate::Putinar{<:Newton},
    domain::BasicSemialgebraicSet,
    p,
)
    return with_fixed_basis(
        domain,
        p,
        certificate.maxdegree,
        certificate.multipliers_certificate.newton,
    )
end

function preorder_indices(
    ::Putinar,
    domain::Union{WithVariables,WithFixedBases},
)
    return map(PreorderIndex, eachindex(domain.inner.p))
end

multiplier_maxdegree(maxdegree, q) = maxdegree - MP.maxdegree(q)
function multiplier_basis(
    certificate::Putinar{<:MaxDegree},
    index::PreorderIndex,
    domain::WithVariables,
)
    q = domain.inner.p[index.value]
    return maxdegree_gram_basis(
        certificate.multipliers_certificate.gram_basis,
        MP.variables(domain),
        # FIXME we ignore `certificate.multipliers_certificate.maxdegree` as it should be the same anyway
        multiplier_maxdegree(certificate.maxdegree, q),
    )
end
function multiplier_basis(
    ::Putinar{<:Newton},
    index::PreorderIndex,
    domain::WithFixedBases,
)
    return domain.bases[index.value]
end
function multiplier_basis_type(::Type{<:Putinar{MC}}, ::Type) where {MC}
    return MA.promote_operation(gram_basis, MC)
end

function generator(
    ::Putinar,
    index::PreorderIndex,
    domain::Union{WithVariables,WithFixedBases},
)
    return domain.inner.p[index.value]
end

ideal_certificate(certificate::Putinar) = certificate.ideal_certificate
ideal_certificate(::Type{<:Putinar{MC,IC}}) where {MC,IC} = IC

"""
    gram_basis(
        certificate::AbstractPreorderCertificate,
        preprocessed,
        poly,
    )

Return the gram basis of `σ_0`, the ideal multiplier in the Putinar-style
decomposition `p = σ_0 + Σ g_i σ_i`.

Whenever the certificate already pre-computes `σ_0`'s basis as a side-effect
of [`preprocessed_domain`](@ref) (e.g. `Putinar{<:Newton}`, which gets it
from [`half_newton_polytope`](@ref)), the preprocessed value is returned
directly and `poly` is ignored. Otherwise (`Putinar{<:MaxDegree}` and
similar), the fallback computes the basis by calling
[`gram_basis`](@ref) on the inner [`ideal_certificate`](@ref) with `poly`.
"""
function gram_basis(::Putinar{<:Newton}, preprocessed::WithFixedBases, _poly)
    return preprocessed.ideal_basis
end

function gram_basis(certificate::Putinar, preprocessed::WithVariables, poly)
    return gram_basis(
        ideal_certificate(certificate),
        with_variables(poly, preprocessed),
    )
end

function SumOfSquares.matrix_cone_type(::Type{<:Putinar{MC}}) where {MC}
    return SumOfSquares.matrix_cone_type(MC)
end
