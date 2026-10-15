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

In addition to the scalar constraints `g_i(x) ≥ 0` and `h_j(x) = 0` of
[`PolyJuMP.AbstractRelaxationOptimizer`](@ref), *polynomial matrix
inequalities* `G_k(x) ⪰ 0` are supported, where the entries of `G_k` are
polynomials in `x`. Their multiplier in the Lagrangian is a Sum-of-Squares
*matrix* `Σ_k` of the same size, the term `Σ_i σ_i * g_i` of the certificate
being completed by `Σ_k ⟨G_k, Σ_k⟩` [HL06]. Such a constraint is added by
constraining a matrix of polynomial expressions to `PSDCone()`; a
[`SOSCone`](@ref) constraint whose coefficients are affine in `x` is also
reformulated into one by the Sum-of-Squares bridges, its Gram matrix becoming
`G_k`.

[HL06] Henrion, Didier and Lasserre, Jean-Bernard.
*Convergent relaxations of polynomial matrix inequalities and static output
feedback.*
IEEE Transactions on Automatic Control 51.2 (2006): 192-202.
"""
mutable struct Optimizer{T} <: PolyJuMP.AbstractRelaxationOptimizer{T}
    model::PolyJuMP.Model{T}
    multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
    matrix_inequalities::Vector{Matrix{PolyJuMP.PolyType{T}}}
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
        Matrix{PolyJuMP.PolyType{T}}[],
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

# The matrix inequalities are not stored in `model.model` so their violation
# is not taken into account by `PolyJuMP.Solution`.
function _add_matrix_violation(
    model::Optimizer{T},
    sol::PolyJuMP.Solution{T},
) where {T}
    if isempty(model.matrix_inequalities)
        return sol
    end
    x = MP.variables(model.model)
    ε = sol.max_constraint_violation
    for G in model.matrix_inequalities
        M = map(g -> convert(T, g(x => sol.values)), G)
        ε = max(ε, -minimum(eigvals(Symmetric(M))))
    end
    tol = model.feasibility_tolerance
    status = if ε < tol
        MOI.FEASIBLE_POINT
    elseif ε < 100tol
        MOI.NEARLY_FEASIBLE_POINT
    else
        MOI.INFEASIBLE_POINT
    end
    return PolyJuMP.Solution{T}(sol.values, sol.objective_value, ε, status)
end

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
        sol =
            PolyJuMP.Solution(values, model.model, model.feasibility_tolerance)
        push!(solutions, _add_matrix_violation(model, sol))
    end
    return solutions
end

# Matrix inequalities
#
# `PolyJuMP.Model` only stores scalar (in)equalities in its
# `SemialgebraicSets` set so the matrix inequalities are stored here instead
# and `PolyJuMP._optimize!` is specialized below to add their multiplier to
# the Lagrangian.

const _TriangleCone = Union{
    MOI.PositiveSemidefiniteConeTriangle,
    PositiveSemidefinite2x2ConeTriangle,
}

const _MatrixCone = Union{_TriangleCone,MOI.PositiveSemidefiniteConeSquare}

function MOI.is_empty(model::Optimizer)
    return MOI.is_empty(model.model) && isempty(model.matrix_inequalities)
end

function MOI.empty!(model::Optimizer)
    MOI.empty!(model.model)
    empty!(model.multiplier_maxdegree)
    empty!(model.matrix_inequalities)
    PolyJuMP._invalidate!(model)
    return
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{
        <:Union{
            MOI.VectorOfVariables,
            MOI.VectorAffineFunction{T},
            MOI.VectorQuadraticFunction{T},
            MOI.VectorNonlinearFunction,
        },
    },
    ::Type{<:_MatrixCone},
) where {T}
    return true
end

# `PolyJuMP.AbstractPolynomialOptimizer` forwards the model attributes to the
# `PolyJuMP.Model` it stores but the bridging costs must be answered by this
# optimizer as it supports more constraints than what that model stores.
function MOI.get(
    model::Optimizer,
    ::MOI.VariableBridgingCost{S},
) where {S<:MOI.AbstractVectorSet}
    return MOI.supports_add_constrained_variables(model, S) ? 0.0 : Inf
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintBridgingCost{F,S},
) where {F<:MOI.AbstractVectorFunction,S<:MOI.AbstractVectorSet}
    return MOI.supports_constraint(model, F, S) ? 0.0 : Inf
end

_scalarize(func::MOI.VectorNonlinearFunction) = func.rows
_scalarize(func::MOI.AbstractVectorFunction) = MOI.Utilities.scalarize(func)

function _entry(model::Optimizer{T}, func) where {T}
    return PolyJuMP._polynomial(
        model.model.variables,
        PolyJuMP._to_polynomial(func, T),
    )
end

# The entries of a triangle cone are given column by column so
# `(1, 1), (1, 2), (2, 2), (1, 3), ...`
function _matrix(model::Optimizer{T}, func, set::_TriangleCone) where {T}
    scalars = _scalarize(func)
    n = MOI.side_dimension(set)
    G = Matrix{PolyJuMP.PolyType{T}}(undef, n, n)
    k = 0
    for j in 1:n, i in 1:j
        k += 1
        G[i, j] = G[j, i] = _entry(model, scalars[k])
    end
    return G
end

# `MOI.PositiveSemidefiniteConeSquare` gives the `n^2` entries column by
# column and the symmetry of the matrix is part of the set so the entries
# that are not syntactically equal give an equality constraint.
function _matrix(
    model::Optimizer{T},
    func,
    set::MOI.PositiveSemidefiniteConeSquare,
) where {T}
    scalars = _scalarize(func)
    n = MOI.side_dimension(set)
    G = Matrix{PolyJuMP.PolyType{T}}(undef, n, n)
    for j in 1:n, i in 1:n
        G[i, j] = _entry(model, scalars[(j-1)*n+i])
    end
    for j in 1:n, i in 1:(j-1)
        if !iszero(G[i, j] - G[j, i])
            model.model.set =
                model.model.set ∩
                SemialgebraicSets.equality(G[i, j] - G[j, i], zero(T))
            G[i, j] = G[j, i] = (G[i, j] + G[j, i]) / 2
        end
    end
    return G
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.AbstractVectorFunction,
    set::_MatrixCone,
)
    push!(model.matrix_inequalities, _matrix(model, func, set))
    PolyJuMP._invalidate!(model)
    i = length(model.matrix_inequalities)
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(i)
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,<:_MatrixCone},
)
    return ci.value in eachindex(model.matrix_inequalities)
end

function MOI.supports(
    ::Optimizer,
    ::PolyJuMP.MultiplierMaxdegree,
    ::Type{<:MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,<:_MatrixCone}},
)
    return true
end

function MOI.set(
    model::Optimizer,
    ::PolyJuMP.MultiplierMaxdegree,
    ci::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,<:_MatrixCone},
    degree::Integer,
)
    MOI.throw_if_not_valid(model, ci)
    model.multiplier_maxdegree[ci] = degree
    PolyJuMP._invalidate!(model)
    return
end

function MOI.get(
    model::Optimizer,
    ::PolyJuMP.MultiplierMaxdegree,
    ci::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,<:_MatrixCone},
)
    return get(model.multiplier_maxdegree, ci, nothing)
end

_maxdegree(G::Matrix) = maximum(MP.maxdegree, G)

function _matrix_multiplier(relaxation, x, n, degree)
    Σ = Matrix{MP.polynomial_type(PolyJuMP.VarType, JuMP.AffExpr)}(undef, n, n)
    for j in 1:n, i in 1:j
        Σ[i, j] = Σ[j, i] = PolyJuMP._multiplier(relaxation, x, degree)
    end
    return Σ
end

# This is a copy of `PolyJuMP._optimize!` with the addition of the multipliers
# of the matrix inequalities; it should be moved to `PolyJuMP` once
# `PolyJuMP.Model` can store matrix inequalities.
function PolyJuMP._optimize!(model::Optimizer{T}) where {T}
    pop = model.model
    x = MP.variables(pop)
    if pop.objective_sense == MOI.FEASIBILITY_SENSE ||
       isnothing(pop.objective_function)
        f = zero(PolyJuMP.PolyType{T})
    else
        f = pop.objective_function
    end
    eqs = PolyJuMP._equalities(pop.set)
    ineqs = PolyJuMP._inequalities(pop.set)
    mats = model.matrix_inequalities
    maxdeg = max(
        MP.maxdegree(f),
        maximum(MP.maxdegree, eqs; init = 0),
        maximum(MP.maxdegree, ineqs; init = 0),
        maximum(_maxdegree, mats; init = 0),
    )
    maxdeg += isodd(maxdeg)
    eq_degrees = Dict{Int,Int}()
    ineq_degrees = Dict{Int,Int}()
    mat_degrees = Dict{Int,Int}()
    for (ci, d) in model.multiplier_maxdegree
        if ci isa MOI.ConstraintIndex{<:Any,<:_MatrixCone}
            mat_degrees[ci.value] = d
        elseif ci isa MOI.ConstraintIndex{<:Any,<:MOI.EqualTo}
            eq_degrees[ci.value] = d
        else
            ineq_degrees[ci.value] = d
        end
    end
    relaxation = JuMP.GenericModel{T}(model.solver)
    cone = PolyJuMP.nonnegativity_cone(model)
    t = JuMP.@variable(relaxation, base_name = "t")
    lagrangian = f - t
    if pop.objective_sense == MOI.MAX_SENSE
        lagrangian = -lagrangian
    end
    for (i, h) in enumerate(eqs)
        d = get(eq_degrees, i, maxdeg - MP.maxdegree(h))
        lagrangian -= PolyJuMP._multiplier(relaxation, x, d) * h
    end
    for (i, g) in enumerate(ineqs)
        d = get(ineq_degrees, i, maxdeg - MP.maxdegree(g))
        σ = PolyJuMP._multiplier(relaxation, x, d)
        JuMP.@constraint(relaxation, σ in cone)
        lagrangian -= σ * g
    end
    for (i, G) in enumerate(mats)
        d = get(mat_degrees, i, maxdeg - _maxdegree(G))
        Σ = _matrix_multiplier(relaxation, x, size(G, 1), d)
        JuMP.@constraint(relaxation, Σ in SOSMatrixCone())
        lagrangian -= sum(G .* Σ)
    end
    cref = JuMP.@constraint(relaxation, lagrangian in cone)
    sense = pop.objective_sense == MOI.MAX_SENSE ? MOI.MIN_SENSE : MOI.MAX_SENSE
    JuMP.set_objective(relaxation, sense, t)
    JuMP.optimize!(relaxation)
    model.relaxation = relaxation
    model.solutions =
        PolyJuMP.recover_solutions(model, relaxation, cref, lagrangian)
    PolyJuMP.postprocess!(model.solutions, pop, nothing)
    return
end
