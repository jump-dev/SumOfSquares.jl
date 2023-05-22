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
