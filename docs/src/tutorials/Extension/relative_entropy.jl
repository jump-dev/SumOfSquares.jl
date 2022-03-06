# # Relative Entropy

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extension/relative_entropy.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extension/relative_entropy.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
using DynamicPolynomials
@polyvar x y
motzkin = x^4*y^2 + x^2*y^4 + 1 - 3x^2*y^2

module RelativeEntropy

import MutableArithmetics
const MA = MutableArithmetics
import MultivariatePolynomials
const MP = MultivariatePolynomials
import MathOptInterface
const MOI = MathOptInterface
import JuMP
import PolyJuMP

# **S**ums of **A**M/**G**M **E**xponential
struct SAGECone <: MOI.AbstractVectorSet
    α::Matrix{Int}
end
MOI.dimension(set::SAGECone) = size(set.α, 1)
Base.copy(set::SAGECone) = set

struct SAGESet <: PolyJuMP.PolynomialSet end
JuMP.reshape_set(::SAGECone, ::PolyJuMP.PolynomialShape) = SAGESet()
function _exponents_matrix(monos)
    α = Matrix{Int}(undef, length(monos), MP.nvariables(monos))
    for (i, mono) in enumerate(monos)
        α[i, :] = MP.exponents(mono)
    end
    return α
end
JuMP.moi_set(set::SAGESet, monos) = SAGECone(_exponents_matrix(monos))

struct AGECone <: MOI.AbstractVectorSet
    α::Matrix{Int}
    k::Int
end
MOI.dimension(set::AGECone) = size(set.α, 1)
Base.copy(set::AGECone) = set

struct AGESet{MT <: MP.AbstractMonomial} <: PolyJuMP.PolynomialSet
    monomial::MT
end
function JuMP.reshape_set(set::AGECone, shape::PolyJuMP.PolynomialShape)
    return AGESet(shape.monomials[set.k])
end
function JuMP.moi_set(set::AGESet, monos)
    k = findfirst(isequal(set.monomial), monos)
    return AGECone(_exponents_matrix(monos), k)
end

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SAGESet)
    return
end

function JuMP.build_constraint(_error::Function, p, set::Union{SAGESet,AGESet}; kws...)
    coefs = PolyJuMP.non_constant_coefficients(p)
    monos = MP.monomials(p)
    cone = JuMP.moi_set(set, monos)
    shape = PolyJuMP.PolynomialShape(monos)
    return PolyJuMP.bridgeable(
        JuMP.VectorConstraint(coefs, cone, shape),
        JuMP.moi_function_type(typeof(coefs)),
        typeof(cone),
    )
end

struct SAGEBridge{T,F,G} <: MOI.Bridges.Constraint.AbstractBridge
    ν::Matrix{MOI.VariableIndex}
    age_constraints::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables,AGECone}}
    equality_constraints::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SAGEBridge{T,F,G}},
    model,
    func::G,
    set::SAGECone,
) where {T,F,G}
    m = size(set.α, 1)
    ν = Vector{Vector{MOI.VariableIndex}}(undef, m)
    age_constraints = Vector{MOI.ConstraintIndex{MOI.VectorOfVariables,AGECone}}(undef, m)
    for k in 1:m
        ν[i], age_constraints[i] = MOI.add_constrained_variables(model, AGECone(set.α, k))
    end
    scalars = MOI.Utilities.eachscalar(func)
    n = size(set.α, 2)
    equality_constraints = map(1:n) do i
        f = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(one(T), ν[:, i]), zero(T))
        MOI.add_constraint(model, MA.sub!(f, scalars[i]), MOI.EqualTo(zero(T)))
    end
    return SAGEBridge{T,F,G}(ν, age_constraints, equality_constraints)
end

function MOI.supports_constraint(
    ::Type{<:SAGEBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SAGECone},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:SAGEBridge})
    return Tuple{Type}[AGECone]
end

function MOI.Bridges.added_constraint_types(::Type{<:SAGEBridge{T,F}}) where {T,F}
    return Tuple{Type,Type}[(F, MOI.EqualTo{T})]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SAGEBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SAGECone},
) where {T}
    S = MOI.Utilities.scalar_type(G)
    F = MOI.Utilities.promote_operation(-, T, MOI.ScalarAffineFunction{Float64}, S)
    return AGEBridge{T,F,G}
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SAGECone},
)
    return [SAGEBridge]
end

struct AGEBridge{T,F,G,H} <: MOI.Bridges.Constraint.AbstractBridge
    k::Int
    equality_constraints::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
    relative_entropy_constraint::MOI.ConstraintIndex{G,MOI.RelativeEntropyCone}
end # See https://arxiv.org/pdf/1907.00814.pdf

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{AGEBridge{T,F,G,H}},
    model,
    func::G,
    set::AGECone,
) where {T,F,G,H}
    m = size(set.α, 1)
    ν = MOI.add_variables(model, m - 1)
    sumν = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(one(T), ν), zero(T))
    ceq = map(1:size(set.α, 2)) do var
        f = -sumν
        j = 0
        for i in 1:m
            if i != set.k
                j += 1
                MA.add_mul!(f, convert(T, set.α[i, var]), MOI.SingleVariable(ν[j]))
            end
        end
        MOI.Utilities.normalize_and_add_constraint(model, f, MOI.EqualTo(zero(T)))
    end
    scalars = MOI.Utilities.eachscalar(func)
    f = MOI.Utilities.operate(
        vcat,
        T,
        scalars[set.k] + sumν,
        scalars[setdiff(1:m, set.k)],
        MOI.VectorOfVariables(ν)
    )
    relative_entropy_constraint = MOI.add_constraint(model, f, MOI.RelativeEntropyCone(2m - 1))
    return AGEBridge{T,F,G,H}(set.k, ceq, relative_entropy_constraint)
end

function MOI.supports_constraint(
    ::Type{<:AGEBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:AGECone},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:AGEBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(::Type{<:AGEBridge{T,F,G}}) where {T,F,G}
    return [(F, MOI.EqualTo{T}), (G, MOI.RelativeEntropyCone)]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:AGEBridge{T}},
    H::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:AGECone},
) where {T}
    S = MOI.Utilities.scalar_type(H)
    F = MOI.Utilities.promote_operation(+, T, S, MOI.ScalarAffineFunction{Float64})
    G = MOI.Utilities.promote_operation(vcat, T, S, F)
    return AGEBridge{T,F,G,H}
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:AGECone},
)
    return [AGEBridge]
end

end

using PolyJuMP
model = Model()
PolyJuMP.setpolymodule!(model, RelativeEntropy)
@constraint(model, motzkin >= 0)

# As we can see above, even if the constraint is internally stored as a
# `MOI.VectorAffineFunction`-in-`SAGECone`, the constraint is reshaped
# into a `Polynomial`-in-`SOSSet` for printing thanks to the `reshape_set` we
# have implemented.

model = Model()
@constraint(model, motzkin in RelativeEntropy.AGESet(x^2*y^2))

import ECOS
set_optimizer(model, ECOS.Optimizer)
optimize!(model)

@show termination_status(model)
