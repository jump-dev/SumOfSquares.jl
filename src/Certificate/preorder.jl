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

struct WithVariables{S,V}
    inner::S
    variables::V
end

function MP.variables(v::WithVariables)
    return v.variables
end
SA.basis(v::WithVariables) = SA.basis(v.inner)
_algebra_element(v::WithVariables) = v.inner
_algebra_element(a::SA.AlgebraElement) = a

struct WithFixedBases{S,B}
    inner::S
    bases::Vector{B}
end

_merge_sorted(a::Vector, ::Tuple{}) = a
function _merge_sorted(a::Vector, b::Vector)
    vars = sort!(vcat(a, b), rev = true)
    unique!(vars)
    return vars
end
_merge_sorted(a::Tuple{}, ::Tuple{}) = a
_merge_sorted(a::Tuple, ::Tuple{}) = a
_merge_sorted(::Tuple{}, b::Tuple) = b
function _merge_sorted(a::Tuple, b::Tuple)
    v = first(a)
    w = first(b)
    if v == w
        return (v, _merge_sorted(Base.tail(a), Base.tail(b))...)
    elseif v > w
        return (v, _merge_sorted(Base.tail(a), b)...)
    else
        return (w, _merge_sorted(a, Base.tail(b))...)
    end
end

_vars(::SemialgebraicSets.FullSpace) = tuple()
function _vars(x::SA.AlgebraElement)
    if SA.basis(x) isa SA.ImplicitBasis
        return MP.variables(SA.coeffs(x))
    else
        return MP.variables(SA.basis(x))
    end
end
_vars(x) = MP.variables(x)

function with_variables(inner, outer)
    return WithVariables(inner, _merge_sorted(_vars(inner), _vars(outer)))
end

function with_fixed_basis(
    domain,
    p,
    maxdegree,
    newton::AbstractNewtonPolytopeApproximation,
)
    v = with_variables(domain, p)
    return WithFixedBases(
        v.inner,
        half_newton_polytope(p, SemialgebraicSets.inequalities(domain), v.variables, maxdegree, newton),
    )
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
        certificate.multipliers_certificate.basis,
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
function multiplier_basis_type(::Type{<:Putinar{MC}}) where {MC}
    return gram_basis_type(MC)
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

function SumOfSquares.matrix_cone_type(::Type{<:Putinar{MC}}) where {MC}
    return SumOfSquares.matrix_cone_type(MC)
end
