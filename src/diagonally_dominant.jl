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
