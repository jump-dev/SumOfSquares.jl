"""
    struct EmptyCone <: MOI.AbstractSymmetricMatrixSetTriangle end

Cone of positive semidefinite matrices of 0 rows and 0 columns. It is reformulated into no
constraints which is useful as some solvers do not support positive semidefinite
cones of dimension zero.
"""
struct EmptyCone <: MOI.AbstractSymmetricMatrixSetTriangle end
MOI.side_dimension(::EmptyCone) = 0

"""
    struct PositiveSemidefinite2x2ConeTriangle <: MOI.AbstractSymmetricMatrixSetTriangle end

Cone of positive semidefinite matrices of 2 rows and 2 columns. It is reformulated into an
`MOI.RotatedSecondOrderCone` constraint, see for instance [Fawzi2018; Equation (1)](@cite).
"""
struct PositiveSemidefinite2x2ConeTriangle <:
       MOI.AbstractSymmetricMatrixSetTriangle end
MOI.side_dimension(::PositiveSemidefinite2x2ConeTriangle) = 2

# isbits types, nothing to copy
function Base.copy(set::Union{EmptyCone,PositiveSemidefinite2x2ConeTriangle})
    return set
end

function matrix_cone(
    S::Type{<:MOI.AbstractSymmetricMatrixSetTriangle},
    side_dimension,
)
    if iszero(side_dimension)
        # Some solvers such as Mosek does not support 0-dimensional PSD cone
        return EmptyCone()
    elseif isone(side_dimension)
        # PSD constraints on 1x1 matrices are equivalent to the nonnegativity
        # of the only entry.
        return MOI.Nonnegatives(1)
    elseif side_dimension == 2
        # PSD constraints on 2x2 matrices are SOC representable.
        return PositiveSemidefinite2x2ConeTriangle()
    else
        # PSD constraints on nxn matrices with n ≥ 3 is not SOC representable,
        # see [F18].
        #
        # [F18] Fawzi, Hamza
        # On representing the positive semidefinite cone using the second-order
        # cone.
        # Mathematical Programming (2018): 1-10.
        return S(side_dimension)
    end
end

"""
    struct DiagonallyDominantConeTriangle <: MOI.AbstractSymmetricMatrixSetTriangle
        side_dimension::Int
    end

See [Ahmadi2017; Definition 4](@cite) for a precise definition of the last two items.
"""
struct DiagonallyDominantConeTriangle <: MOI.AbstractSymmetricMatrixSetTriangle
    side_dimension::Int
end

function matrix_cone(S::Type{DiagonallyDominantConeTriangle}, side_dimension)
    if iszero(side_dimension)
        return EmptyCone()
    elseif isone(side_dimension)
        return MOI.Nonnegatives(1)
    else
        # With `side_dimension` = 2, we want to avoid using SOC and only use LP
        return S(side_dimension)
    end
end

"""
    struct ScaledDiagonallyDominantConeTriangle <: MOI.AbstractSymmetricMatrixSetTriangle
        side_dimension::Int
    end

See [Ahmadi2017; Definition 4](@cite) for a precise definition of the last two items.
"""
struct ScaledDiagonallyDominantConeTriangle <:
       MOI.AbstractSymmetricMatrixSetTriangle
    side_dimension::Int
end

function matrix_cone(
    S::Type{ScaledDiagonallyDominantConeTriangle},
    side_dimension,
)
    if iszero(side_dimension)
        return EmptyCone()
    elseif isone(side_dimension)
        return MOI.Nonnegatives(side_dimension)
    elseif side_dimension == 2
        return PositiveSemidefinite2x2ConeTriangle()
    else
        return S(side_dimension)
    end
end

# isbits types, nothing to copy
function Base.copy(
    set::Union{
        DiagonallyDominantConeTriangle,
        ScaledDiagonallyDominantConeTriangle,
    },
)
    return set
end

"""
    struct WeightedSOSCone{B,G,W}
        basis::B
        gram_bases::Vector{G}
        weights::Vector{W}
    end

The weighted sum-of-squares cone is the set of vectors of coefficients `a` in `basis`
that are the sum of `weights[i]` multiplied by a gram matrix with basis `gram_bases[i]`.

See [Papp2017; Section 1.1](@cite) and [Kapelevich2023; Section 1](@cite).
"""
struct WeightedSOSCone{
    M,
    B<:AbstractPolynomialBasis,
    G<:AbstractPolynomialBasis,
    W<:MP.AbstractPolynomialLike,
} <: MOI.AbstractVectorSet
    basis::B
    gram_bases::Vector{G}
    weights::Vector{W}
end
MOI.dimension(set::WeightedSOSCone) = length(set.basis)
Base.copy(set::WeightedSOSCone) = set

struct SOSPolynomialSet{
    DT<:AbstractSemialgebraicSet,
    B<:AbstractPolynomialBasis,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
    CT<:Certificate.AbstractCertificate,
} <: MOI.AbstractVectorSet
    domain::DT
    basis::B
    certificate::CT
end
MOI.dimension(set::SOSPolynomialSet) = length(set.monomials)
Base.copy(set::SOSPolynomialSet) = set
