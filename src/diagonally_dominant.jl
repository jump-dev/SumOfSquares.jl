function matrix_cone(S::Type{<:Union{MatrixConeTriangle,
                                     MOI.PositiveSemidefiniteConeTriangle}},
                     side_dimension)
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
    struct DiagonallyDominantConeTriangle <: MatrixConeTriangle
        side_dimension::Int
    end

See Definition 4 of [AM17] for a precise definition of the last two items.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.
"""
struct DiagonallyDominantConeTriangle <: MatrixConeTriangle
    side_dimension::Int
end

function matrix_cone(S::Type{DiagonallyDominantConeTriangle},
                     side_dimension)
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
    struct ScaledDiagonallyDominantConeTriangle <: MatrixConeTriangle
        side_dimension::Int
    end

See Definition 4 of [AM17] for a precise definition of the last two items.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.
"""
struct ScaledDiagonallyDominantConeTriangle <: MatrixConeTriangle
    side_dimension::Int
end

function matrix_cone(S::Type{ScaledDiagonallyDominantConeTriangle},
                     side_dimension)
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

function side_dimension(set::Union{MatrixConeTriangle,
                                   MOI.PositiveSemidefiniteConeTriangle})
    return set.side_dimension
end

# isbits types, nothing to copy
function Base.copy(set::MatrixConeTriangle)
    return set
end

# TODO make PSDConeTriangle inherit from MatrixConeTriangle to remove the need
#      for this
function MOI.dimension(set::MatrixConeTriangle)
    return div(side_dimension(set) * (side_dimension(set) + 1), 2)
end

function MOIU.set_dot(x::Vector, y::Vector, set::MatrixConeTriangle)
    return MOIU.triangle_dot(x, y, side_dimension(set), 0)
end

function MOIU.dot_coefficients(a::Vector, set::MatrixConeTriangle)
    b = copy(a)
    MOIU.triangle_coefficients!(b, side_dimension(set), 0)
    return b
end
