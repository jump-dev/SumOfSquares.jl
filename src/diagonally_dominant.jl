abstract type MatrixConeTriangle <: MOI.AbstractVectorSet end

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

# isbits types, nothing to copy
function Base.copy(set::MatrixConeTriangle)
    return set
end

# TODO make PSDConeTriangle inherit from MatrixConeTriangle to remove the need
#      for this
function MOI.dimension(set::MatrixConeTriangle)
    return div(set.side_dimension * (set.side_dimension + 1), 2)
end

function MOIU.set_dot(x::Vector, y::Vector, set::MatrixConeTriangle)
    return MOIU.triangle_dot(x, y, set.side_dimension, 0)
end

function MOIU.dot_coefficients(a::Vector, set::MatrixConeTriangle)
    b = copy(a)
    MOIU.triangle_coefficients!(b, set.side_dimension, 0)
    return b
end
