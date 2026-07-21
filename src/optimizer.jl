"""
    Optimizer{T}(solver)

Optimizer computing a bound on the objective value of a polynomial
optimization problem using its Lasserre / Sum-of-Squares relaxation, see
[`PolyJuMP.AbstractRelaxationOptimizer`](@ref) with the cone [`SOSCone`](@ref)
as certificate of nonnegativity.
The level of the hierarchy is chosen for each constraint with the
`PolyJuMP.MultiplierMaxdegree` constraint attribute.
The relaxation is solved with `solver` and the bound can be queried with
`MOI.ObjectiveBound()`.
"""
mutable struct Optimizer{T} <: PolyJuMP.AbstractRelaxationOptimizer{T}
    model::PolyJuMP.Model{T}
    multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
    solver::Any
    relaxation::Union{Nothing,JuMP.GenericModel{T}}
    solutions::Vector{PolyJuMP.Solution{T}}
    feasibility_tolerance::T
    solve_time::Float64
end

function Optimizer{T}(
    solver;
    feasibility_tolerance = sqrt(Base.rtoldefault(T)),
) where {T}
    return Optimizer{T}(
        PolyJuMP.Model{T}(),
        Dict{MOI.ConstraintIndex,Int}(),
        solver,
        nothing,
        PolyJuMP.Solution{T}[],
        feasibility_tolerance,
        NaN,
    )
end

Optimizer(solver; kws...) = Optimizer{Float64}(solver; kws...)

MOI.get(::Optimizer, ::MOI.SolverName) = "SumOfSquares"

PolyJuMP.nonnegativity_cone(::Optimizer) = SOSCone()

# Solutions could be recovered from the moment matrix given by the dual of
# the SOS constraint using `MultivariateMoments.atomic_measure`
function PolyJuMP.recover_solutions(
    ::Optimizer{T},
    ::JuMP.GenericModel{T},
    ::JuMP.ConstraintRef,
    lagrangian,
) where {T}
    return PolyJuMP.Solution{T}[]
end
