struct EmptyCone <: MOI.AbstractSymmetricMatrixSetTriangle end
MOI.side_dimension(::EmptyCone) = 0

"""
    struct PositiveSemidefinite2x2ConeTriangle <: MOI.AbstractSymmetricMatrixSetTriangle end

Cone of positive semidefinite matrices of 2 rows and 2 columns.
"""
struct PositiveSemidefinite2x2ConeTriangle <:
       MOI.AbstractSymmetricMatrixSetTriangle end
MOI.side_dimension(::PositiveSemidefinite2x2ConeTriangle) = 2

# isbits types, nothing to copy
function Base.copy(set::Union{EmptyCone,PositiveSemidefinite2x2ConeTriangle})
    return set
end
