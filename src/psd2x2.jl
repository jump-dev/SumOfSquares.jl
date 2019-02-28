abstract type MatrixConeTriangle <: MOI.AbstractVectorSet end

struct EmptyCone <: MatrixConeTriangle end
side_dimension(::EmptyCone) = 0

struct PositiveSemidefinite2x2ConeTriangle <: MatrixConeTriangle end
side_dimension(::PositiveSemidefinite2x2ConeTriangle) = 2
