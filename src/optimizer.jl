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

"""
    PolyJuMP.recover_solutions(model::Optimizer, relaxation, cref, lagrangian)

Recover candidate solutions from the dual of the SOS constraint `cref` of the
Lagrangian. This dual is a moment matrix and, if the relaxation is tight and
the moments correspond to an atomic measure, the atoms of this measure are
optimal solutions [HL05]. The extraction of the atoms is implemented by
`MultivariateMoments.atomic_measure`; an empty vector of solutions is
returned when it detects that the moment matrix is not atomic.

[HL05] Henrion, Didier, and Jean-Bernard Lasserre.
"Detecting global optimality and extracting solutions in GloptiPoly."
Positive polynomials in control. Springer (2005): 293-310.
"""
function PolyJuMP.recover_solutions(
    model::Optimizer{T},
    relaxation::JuMP.GenericModel{T},
    cref::JuMP.ConstraintRef,
    lagrangian,
) where {T}
    solutions = PolyJuMP.Solution{T}[]
    if JuMP.termination_status(relaxation) != MOI.OPTIMAL
        return solutions
    end
    ν = MultivariateMoments.moment_matrix(cref)
    measure = MultivariateMoments.atomic_measure(ν, sqrt(Base.rtoldefault(T)))
    if isnothing(measure)
        return solutions
    end
    x = MP.variables(model.model)
    for atom in measure.atoms
        values = zeros(T, length(x))
        for (j, var) in enumerate(measure.variables)
            values[findfirst(isequal(var), x)] = atom.center[j]
        end
        push!(
            solutions,
            PolyJuMP.Solution(values, model.model, model.feasibility_tolerance),
        )
    end
    return solutions
end
