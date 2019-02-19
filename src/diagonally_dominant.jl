abstract type MatrixConeTriangle <: MOI.AbstractVectorSet end

function matrix_cone(S::Type{<:Union{MatrixConeTriangle,
                                     MOI.PositiveSemidefiniteConeTriangle}},
                     side_dimension)
    return S(side_dimension)
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
