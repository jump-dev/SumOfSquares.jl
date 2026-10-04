# Heuristics to obtain feasible points from the moments computed by the
# Sum-of-Squares program, see the "Rounding" section of the documentation.

"""
    abstract type AbstractRounding end

Strategy used by [`round_solution`](@ref) to generate candidate points from
the moments of a measure. The candidates are then projected onto the feasible
set with [`heuristic_projection`](@ref).
"""
abstract type AbstractRounding end

"""
    FirstMomentRounding()

The only candidate is the vector of first-order moments
``(\\mathbb{E}[x_1], \\ldots, \\mathbb{E}[x_n])``.
If the measure is a Dirac measure at a minimizer, this is the minimizer.
"""
struct FirstMomentRounding <: AbstractRounding end

"""
    GaussianRounding(; num_samples::Int = 100, rng = Random.default_rng())

The candidates are the vector ``\\mu`` of first-order moments followed by
`num_samples` samples of the Gaussian distribution ``\\mathcal{N}(\\mu, \\Sigma)``
where ``\\Sigma = \\mathbb{E}[x x^\\top] - \\mu\\mu^\\top`` is the covariance
matrix of the measure.
This Gaussian has the same moments of order up to 2 as the measure.
For the MAX-CUT problem, this generalizes the random hyperplane rounding of
Goemans and Williamson.
"""
struct GaussianRounding{R<:Random.AbstractRNG} <: AbstractRounding
    num_samples::Int
    rng::R
end

function GaussianRounding(; num_samples::Int = 100, rng = Random.default_rng())
    return GaussianRounding(num_samples, rng)
end

_moment_vector(μ::MomentVector) = μ
_moment_vector(ν::MomentMatrix) = moment_vector(ν)

# Returns the moments of order up to 2 normalized by the mass of the measure.
# The entries corresponding to moments that are not available are `NaN`.
function _low_order_moments(μ::MomentVector{S}, vars) where {S}
    T = float(S)
    n = length(vars)
    mass = zero(T)
    first = fill(T(NaN), n)
    second = fill(T(NaN), n, n)
    for m in moments(μ)
        mono = MP.monomial(m.polynomial)
        d = MP.degree(mono)
        if d > 2
            continue
        end
        I = Int[]
        for (i, var) in enumerate(vars)
            for _ in 1:MP.degree(mono, var)
                push!(I, i)
            end
        end
        if length(I) != d
            error("Variable of monomial `$mono` is not in `$vars`.")
        end
        if d == 0
            mass = moment_value(m)
        elseif d == 1
            first[I[1]] = moment_value(m)
        else
            second[I[1], I[2]] = second[I[2], I[1]] = moment_value(m)
        end
    end
    if !(mass > 0)
        error(
            "Cannot round a measure of mass `$mass`, the solver may have failed.",
        )
    end
    if any(isnan, first)
        error("Missing first-order moments in `$μ`.")
    end
    return first / mass, second / mass
end

"""
    rounding_candidates(μ, vars, rounding::AbstractRounding)

Return a vector of candidate points, i.e., vectors of values for the
variables `vars` built from the moments of `μ`.
"""
function rounding_candidates end

function rounding_candidates(μ, vars, ::FirstMomentRounding)
    m, _ = _low_order_moments(_moment_vector(μ), vars)
    return [m]
end

function rounding_candidates(μ, vars, rounding::GaussianRounding)
    m, M = _low_order_moments(_moment_vector(μ), vars)
    if any(isnan, M)
        error("Missing second-order moments in `$μ`.")
    end
    F = eigen(Symmetric(M - m * m'))
    # The moment matrix computed by the solver may be slightly indefinite
    L = F.vectors * Diagonal(sqrt.(max.(F.values, 0)))
    n = length(vars)
    samples = [m + L * randn(rounding.rng, n) for _ in 1:rounding.num_samples]
    return [[m]; samples]
end

_equalities(K::AbstractSemialgebraicSet) = equalities(K)
_inequalities(::AbstractAlgebraicSet) = []
_inequalities(K::BasicSemialgebraicSet) = inequalities(K)

struct _Constraint
    polynomial::Any
    is_equality::Bool
    # Indices of the variables in `polynomial`
    variables::Vector{Int}
end

mutable struct _Projection{V,T}
    variables::Vector{V}
    constraints::Vector{_Constraint}
    tol::T
    # Number of variables we can still fix at their current value
    max_guesses::Int
end

function _Projection(K, vars, tol, max_guesses)
    polys = [_equalities(K); _inequalities(K)]
    constraints = _Constraint[]
    for (i, p) in enumerate(polys)
        idx = map(MP.effective_variables(p)) do var
            j = findfirst(isequal(var), vars)
            if isnothing(j)
                error("Variable `$var` of the set is not in `$vars`.")
            end
            return j
        end
        is_equality = i <= length(_equalities(K))
        push!(constraints, _Constraint(p, is_equality, idx))
    end
    return _Projection(vars, constraints, tol, max_guesses)
end

function _is_satisfied(proj::_Projection, c::_Constraint, value)
    if c.is_equality
        return abs(value) <= proj.tol
    else
        return value >= -proj.tol
    end
end

function _is_satisfied(proj::_Projection, c::_Constraint, x::Vector)
    return _is_satisfied(proj, c, c.polynomial(proj.variables => x))
end

function _num_free(c::_Constraint, fixed)
    return count(i -> !fixed[i], c.variables)
end

# Coefficients `a` such that `p = sum(a[k + 1] * var^k)`.
function _univariate_coefficients(p, var)
    a = zeros(Float64, MP.maxdegree(p, var) + 1)
    for t in MP.terms(p)
        a[MP.degree(MP.monomial(t), var)+1] += MP.coefficient(t)
    end
    return a
end

_evaluate(a::Vector, t) = evalpoly(t, a)

function _derivative(a::Vector)
    return [k * a[k+1] for k in 1:(length(a)-1)]
end

# Real parts of the roots refined by a few Newton steps. Roots that are not
# real are filtered out later as they will not satisfy the constraints.
function _approximate_real_roots(a::Vector{T}) where {T}
    scale = maximum(abs, a; init = zero(T))
    d = findlast(c -> abs(c) > Base.rtoldefault(T) * scale, a)
    if isnothing(d) || d == 1
        return T[]
    end
    C = zeros(T, d - 1, d - 1)
    for i in 2:(d-1)
        C[i, i-1] = one(T)
    end
    for i in 1:(d-1)
        C[i, end] = -a[i] / a[d]
    end
    roots = real.(eigvals(C))
    da = _derivative(a)
    for i in eachindex(roots)
        for _ in 1:3
            t = roots[i]
            q = _evaluate(a, t)
            dq = _evaluate(da, t)
            if iszero(dq)
                break
            end
            # Newton's method is unstable close to multiple roots so we only
            # accept steps that improve the residual.
            new_t = t - q / dq
            if abs(_evaluate(a, new_t)) >= abs(q)
                break
            end
            roots[i] = new_t
        end
    end
    return roots
end

# Project the `i`th variable onto the set of values satisfying the constraints
# that have no free variable except the `i`th one.
# The closest point to `x[i]` is either `x[i]` or a root of one of the
# univariate polynomials.
function _project_univariate(proj::_Projection, x, fixed, i)
    var = proj.variables[i]
    others = [j for j in eachindex(x) if j != i]
    polys = Vector{Float64}[]
    constraints = _Constraint[]
    for c in proj.constraints
        if i in c.variables && _num_free(c, fixed) == 1
            q = MP.subs(c.polynomial, proj.variables[others] => x[others])
            push!(polys, _univariate_coefficients(q, var))
            push!(constraints, c)
        end
    end
    candidates = [x[i]]
    for a in polys
        append!(candidates, _approximate_real_roots(a))
    end
    best = nothing
    for t in candidates
        if all(eachindex(polys)) do k
            return _is_satisfied(proj, constraints[k], _evaluate(polys[k], t))
        end
            if isnothing(best) || abs(t - x[i]) < abs(best - x[i])
                best = t
            end
        end
    end
    if isnothing(best)
        return
    end
    # Avoid `-0.0`
    return best + zero(best)
end

function _project!(proj::_Projection, x::Vector, fixed::BitVector)
    while !all(fixed)
        # Find a variable that is the only free variable of a constraint
        c = findfirst(c -> _num_free(c, fixed) == 1, proj.constraints)
        if isnothing(c)
            break
        end
        i = proj.constraints[c].variables[findfirst(
            i -> !fixed[i],
            proj.constraints[c].variables,
        )]
        t = _project_univariate(proj, x, fixed, i)
        if isnothing(t)
            return
        end
        x[i] = t
        fixed[i] = true
    end
    # Every constraint has either zero or at least two free variables.
    # We pick the constraint with the fewest free variables and fix one of
    # them at its current value. This makes progress towards univariate
    # constraints. If it fails, we try the next one.
    best = nothing
    for c in proj.constraints
        n = _num_free(c, fixed)
        if n >= 2 && (isnothing(best) || n < _num_free(best, fixed))
            best = c
        end
    end
    if isnothing(best)
        # The free variables do not appear in any constraint
        fixed .= true
        return x
    end
    for i in best.variables
        if fixed[i]
            continue
        end
        if proj.max_guesses <= 0
            return
        end
        proj.max_guesses -= 1
        # Since no constraint has only one free variable, no constraint
        # becomes without free variable when `i` is fixed so there is no
        # constraint to check.
        new_fixed = copy(fixed)
        new_fixed[i] = true
        ret = _project!(proj, copy(x), new_fixed)
        if !isnothing(ret)
            return ret
        end
    end
    return
end

"""
    heuristic_projection(
        K::AbstractSemialgebraicSet,
        vars,
        x::AbstractVector;
        tol = 1e-8,
        max_guesses::Int = 10length(vars),
    )

Heuristically finds a point of `K` close to the point `x` of values for the
variables `vars`. Returns `nothing` if no point was found.

The variables are fixed one by one. If a constraint has only one variable that
is not fixed yet, this variable is set to the closest value to its current value
that satisfies all constraints in which it is the only non-fixed variable.
This is cheap as it amounts to finding the real roots of univariate polynomials.
Otherwise, a variable of a constraint with the fewest non-fixed
variables is fixed at its current value, backtracking if this leads to an
infeasible univariate subproblem. At most `max_guesses` such guesses are tried.
Any point that is returned satisfies the equalities `p(x) = 0` with
`abs(p(x)) <= tol` and the inequalities `p(x) >= 0` with `p(x) >= -tol`.

For instance, if `K` is defined by `0 <= x <= 1` and `y <= x^2`, `x` is first
projected onto `[0, 1]` and then `y` is projected onto `(-∞, x^2]`.
"""
function heuristic_projection(
    K::AbstractSemialgebraicSet,
    vars,
    x::AbstractVector;
    tol = 1e-8,
    max_guesses::Int = 10length(vars),
)
    proj = _Projection(K, vars, tol, max_guesses)
    y = _project!(proj, convert(Vector{Float64}, x), falses(length(vars)))
    if isnothing(y) || !all(c -> _is_satisfied(proj, c, y), proj.constraints)
        return
    end
    return y
end

"""
    round_solution(
        μ,
        K::AbstractSemialgebraicSet,
        objective = nothing;
        rounding::AbstractRounding = FirstMomentRounding(),
        kws...,
    )

Heuristically finds a feasible point of `K` from the moments of the measure `μ`,
which is usually obtained with [`moment_matrix`](@ref) after solving a
Sum-of-Squares program.
Candidate points are generated by `rounding` and are projected onto `K` with
[`heuristic_projection`](@ref) to which the keyword arguments `kws` are
passed.
If `objective` is `nothing`, the first feasible point is returned.
Otherwise, the feasible point minimizing `objective` is returned.
Returns `nothing` if no feasible point was found.
The point is a vector of values for the variables `variables(μ)`, sorted in
the same order as the center of the atoms returned by `atomic_measure`.

Unlike `atomic_measure`, this does not require the moment matrix to be
flat. However, the point found is not necessarily a global minimizer, it only
provides an upper bound on the optimal objective value of the polynomial
optimization problem.
"""
function round_solution(
    μ,
    K::AbstractSemialgebraicSet,
    objective = nothing;
    rounding::AbstractRounding = FirstMomentRounding(),
    kws...,
)
    vars = MP.variables(μ)
    best = nothing
    best_value = Inf
    for x in rounding_candidates(μ, vars, rounding)
        y = heuristic_projection(K, vars, x; kws...)
        if isnothing(y)
            continue
        end
        if isnothing(objective)
            return y
        end
        value = objective(vars => y)
        if value < best_value
            best = y
            best_value = value
        end
    end
    return best
end
