function is_commutative(vars)
    return length(vars) < 2 || prod(vars[1:2]) == prod(reverse(vars[1:2]))
end

abstract type AbstractNewtonPolytopeApproximation end

# `maxdegree` and `mindegree` on each variable separately
# and on each group of variable.
struct NewtonDegreeBounds{NPT} <: AbstractNewtonPolytopeApproximation
    variable_groups::NPT
end

# Filters out points ouside the Newton polytope from the
# outer approximation given by `outer_approximation`.
struct NewtonFilter{N<:AbstractNewtonPolytopeApproximation} <:
       AbstractNewtonPolytopeApproximation
    outer_approximation::N
end

struct DegreeBounds{M}
    mindegree::Int
    maxdegree::Int
    variablewise_mindegree::M
    variablewise_maxdegree::M
end

function __chip(cur, i, vars, exps, n, op)
    if n > 0
        exp = min(exps[i], n)
        next = iszero(exp) ? cur : op(cur, vars[i]^exp)
        return __chip(next, i + 1, vars, exps, n - exp, op)
    else
        return cur
    end
end
function _chip(mono, n)
    vars = MP.variables(mono)
    exps = MP.exponents(mono)
    if n < 0
        vars = reverse(vars)
        exps = reverse(exps)
        op(a, b) = b * a
    else
        op = *
    end
    return __chip(MP.constant_monomial(mono), 1, vars, exps, abs(n), op)
end

function _is_hermitian_square(mono)
    d = MP.degree(mono)
    if isodd(d)
        return false
    end
    n = div(d, 2)
    return _chip(mono, n) == _chip(mono, -n)
end

# TODO add to MP
function _map_powers(f, mono)
    exps = map(f, MP.powers(mono))
    return _monomial(MP.variables(mono), exps)
end
function _map_exponents(f, mono)
    exps = map(f, MP.exponents(mono))
    return _monomial(MP.variables(mono), exps)
end

_min_half(d::Integer) = cld(d, 2)
_max_half(d::Integer) = fld(d, 2)
_half(::Nothing) = nothing
function _half(d::DegreeBounds)
    return DegreeBounds(
        _min_half(d.mindegree),
        _max_half(d.maxdegree),
        _map_exponents(_min_half, d.variablewise_mindegree),
        _map_exponents(_max_half, d.variablewise_maxdegree),
    )
end

function min_degree(p::SA.AlgebraElement, v)
    return mapreduce(Base.Fix2(min_degree, v), min, SA.supp(p))
end
function min_degree(p::MB.Polynomial{B}, v) where {B}
    if _is_monomial_basis(B)
        return MP.degree(p.monomial, v)
    else
        error("TODO $B")
    end
end

function _monomial(vars, exps)
    if any(Base.Fix2(isless, 0), exps)
        return
    end
    return prod(vars .^ exps)
end

function max_degree(p::SA.AlgebraElement, v)
    return mapreduce(Base.Fix2(max_degree, v), max, SA.supp(p))
end
max_degree(p::MB.Polynomial, v) = MP.degree(p.monomial, v)

function min_shift(d, shift)
    return max(0, d - shift)
end

function minus_shift(
    deg,
    mono::MP.AbstractMonomial,
    p::SA.AlgebraElement,
    shift,
)
    return _map_powers(mono) do (v, d)
        return shift(d, deg(p, v))
    end
end

function minus_shift(d::DegreeBounds, p::SA.AlgebraElement)
    var_mindegree =
        minus_shift(min_degree, d.variablewise_mindegree, p, min_shift)
    var_maxdegree = minus_shift(max_degree, d.variablewise_maxdegree, p, -)
    if isnothing(var_maxdegree)
        return
    end
    vars = MP.variables(d.variablewise_maxdegree)
    return DegreeBounds(
        min_shift(d.mindegree, _mindegree(p, vars)),
        d.maxdegree - _maxdegree(p, vars),
        var_mindegree,
        var_maxdegree,
    )
end

function _combine_sign(a, b)
    if ismissing(a) || ismissing(b)
        return missing
    else
        return a == b ? a : missing
    end
end

_sign(a::Number) = sign(a)
# Can be for instance a JuMP or MOI function so the sign can be anything
_sign(_) = missing

function deg_sign(deg, p, d)
    sgn = nothing
    for (k, v) in SA.nonzero_pairs(SA.coeffs(p))
        if deg(SA.basis(p)[k]) == d
            s = _sign(v)
            if isnothing(sgn)
                sgn = s
            else
                sgn = _combine_sign(sgn, s)
            end
        end
    end
    return sgn
end

function deg_sign(deg, p)
    d = deg(p)
    return d, deg_sign(deg, p, d)
end

function _interval(a, b)
    if a < b
        return a:b
    else
        return b:a
    end
end

function _gram_shift_min(d_g, d_s, g::SA.AlgebraElement)
    if _is_monomial_basis(typeof(g))
        return d_g + 2minimum(d_s)
    else
        return 0
    end
end

function deg_range(deg, p, gs, gram_deg, truncation)
    d_max = min(deg(p), truncation)
    sign = deg_sign(deg, p, d_max)
    for g in gs
        d_g, sign_g = deg_sign(deg, g)
        d_s = gram_deg(g)
        if isempty(d_s) || _gram_shift_min(d_g, d_s, g) > truncation
            continue
        end
        d = d_g + 2maximum(d_s)
        if d > truncation
            d = truncation
            if mod(d_g, 2) != mod(truncation, 2)
                d -= 1
            end
        end
        # Multiply by `-1` because we move it to lhs
        # p = s_0 + sum s_i g_i -> p - sum s_i g_i = s_0
        sign_g = -sign_g
        if isnothing(sign)
            d_max = d
            sign = sign_g
        elseif d_max <= d
            if d_max == d
                sign = _combine_sign(sign, sign_g)
            else
                sign = sign_g
            end
            d_max = d
        end
    end
    if !isnothing(sign) && (ismissing(sign) || (iseven(d_max) && sign == 1))
        return true, d_max
    else
        return false, d_max - 1
    end
end

#    deg_range(deg, p, gs, gram_deg, range)
#
# Maximum value of `deg(s_0 = p - sum s_i g_i for i in eachindex(gs)) in range` where
# `s_0, s_i` are SOS and `deg(s_i) <= gram_deg(g)`.
# Note that `range` should be in increasing order.
function deg_range(deg, p, gs, gram_deg, range::UnitRange)
    d = maximum(range)
    while d in range
        ok, d = deg_range(deg, p, gs, gram_deg, d)
        if ok
            return d
        end
    end
    return
end

#_mindegree(p::MP.AbstractPolynomialLike) = MP.mindegree(p)
function _mindegree(a::SA.AlgebraElement, vars)
    return minimum(Base.Fix2(_mindegree, vars), SA.supp(a))
end
function _maxdegree(a::SA.AlgebraElement, vars)
    return maximum(Base.Fix2(_maxdegree, vars), SA.supp(a), init = 0)
end
#_mindegree(basis::MB.SubBasis{MB.Monomial}) = MP.mindegree(basis.monomials)
#_maxdegree(basis::MB.SubBasis{MB.Monomial}) = MP.maxdegree(basis.monomials)
function _mindegree(p::MB.Polynomial{B}, vars) where {B}
    if _is_monomial_basis(B)
        _sum_degree(p.monomial, vars)
    else
        error("TODO $B")
    end
end
_maxdegree(p::MB.Polynomial, vars) = _sum_degree(p.monomial, vars)
function _degree(mono, var::MP.AbstractVariable, comm::Bool)
    if comm
        return MP.degree(mono, var)
    else
        vars = MP.variables(mono)
        return mapreduce(
            j -> vars[j] == var ? MP.exponents(mono)[j] : 0,
            +,
            eachindex(vars),
            init = 0,
        )
    end
end
function _sum_degree(mono, vars)
    comm = is_commutative(vars)
    return sum(var -> _degree(mono, var, comm), vars)
end

_is_monomial_basis(::Type{<:MB.AbstractMonomialIndexed}) = false
_is_monomial_basis(::Type{<:Union{MB.Monomial,MB.ScaledMonomial}}) = true
function _is_monomial_basis(
    ::Type{<:SA.AlgebraElement{<:MB.Algebra{BT,B}}},
) where {BT,B}
    return _is_monomial_basis(B)
end

# Minimum degree of a gram basis for a gram matrix `s`
# such that the minimum degree of `s * g` is at least `mindegree`.
function _multiplier_mindegree(mindegree, g::SA.AlgebraElement, vars)
    if _is_monomial_basis(typeof(g))
        return _min_half(min_shift(mindegree, _mindegree(g, vars)))
    else
        # For instance, with `Chebyshev` a square already produces
        # a monomial of degree 0 so let's just give up here
        # The Chebyshev Newton polytope of the product of polynomials
        # with Chebyshev Newton polytope `N1` and `N2` is a subset
        # of `(N1 + N2 - R_+^n) ∩ R_+^n` so we compute
        # `σ(y, cheby(N1) + cheby(N2))` as `σ(max.(y, 0), N1 + N2)`
        return 0
    end
end

# Maximum degree of a gram basis for a gram matrix `s`
# such that the maximum degree of `s * g` is at most `maxdegree`.
# Here, we assume that the degree of `*(a::MB.Polynomial, b::MB.Polynomial)`
# is bounded by the sum of the degree of `a` and `b` for any basis.
function _multiplier_maxdegree(maxdegree, g::SA.AlgebraElement, vars)
    return _max_half(maxdegree - _maxdegree(g, vars))
end

function _multiplier_deg_range(range, g::SA.AlgebraElement, vars)
    return _multiplier_mindegree(minimum(range), g, vars):_multiplier_maxdegree(
        maximum(range),
        g,
        vars,
    )
end

# Cheap approximation of the convex hull as the approximation of:
#
# z such that mindegree < sum(z) < maxdegree
# |\
# |#\ <-------- sum(z) = maxdegree
# |##\
# |\<-\-------- sum(z) = mindegree
# | \##\
# +---------
#
# and:
#
# z such that minmultideg < z < maxmultideg
# | +----+ <--- maxmultidegree
# | |####|
# | |####|
# | +----+
# | ^---------- minmultidegree
# +---------
function putinar_degree_bounds(
    p::SA.AlgebraElement,
    gs::AbstractVector{<:SA.AlgebraElement},
    vars,
    maxdegree,
)
    mindegree = 0
    # TODO homogeneous case
    degrange(g) = _multiplier_deg_range(mindegree:maxdegree, g, vars)
    minus_degrange(g) = (-).(degrange(g))
    # The multiplier will have degree `0:2fld(maxdegree - _maxdegree(g), 2)`
    mindegree = if _is_monomial_basis(typeof(p))
        d = deg_range(
            (-) ∘ Base.Fix2(_mindegree, vars),
            p,
            gs,
            minus_degrange,
            -maxdegree:0,
        )
        isnothing(d) ? d : -d
    else
        0
    end
    if isnothing(mindegree)
        return
    end
    maxdegree =
        deg_range(Base.Fix2(_maxdegree, vars), p, gs, degrange, 0:maxdegree)
    if isnothing(maxdegree)
        return
    end
    vars_mindeg = map(vars) do v
        return if _is_monomial_basis(eltype(gs))
            -deg_range(
                (-) ∘ Base.Fix2(min_degree, v),
                p,
                gs,
                minus_degrange,
                -maxdegree:0,
            )
        else
            0
        end
    end
    if any(isnothing, vars_mindeg)
        return
    end
    vars_maxdeg = map(vars) do v
        return deg_range(Base.Fix2(max_degree, v), p, gs, degrange, 0:maxdegree)
    end
    if any(isnothing, vars_maxdeg)
        return
    end
    @assert all(d -> d >= 0, vars_mindeg)
    @assert all(d -> d >= 0, vars_maxdeg)
    return DegreeBounds(
        mindegree,
        maxdegree,
        _monomial(vars, vars_mindeg),
        _monomial(vars, vars_maxdeg),
    )
end

function multiplier_basis(
    g::SA.AlgebraElement{<:MB.Algebra{BT,B}},
    bounds::DegreeBounds,
) where {BT,B}
    return maxdegree_gram_basis(
        MB.FullBasis{B,MP.monomial_type(typeof(g))}(),
        _half(minus_shift(bounds, g)),
    )
end

# Cartesian product of the newton polytopes of the different parts
function _cartesian_product(bases::Vector{<:MB.SubBasis{B}}, bounds) where {B}
    monos = [b.monomials for b in bases]
    basis = MB.SubBasis{B}(
        vec([prod(monos) for monos in Iterators.product(monos...)]),
    )
    # We know that the degree inequalities are satisfied variable-wise and
    # part-wise but for all variables together so we filter with that
    if isnothing(bounds)
        return MB.empty_basis(typeof(basis))
    else
        return MB.SubBasis{B}(
            filter(Base.Fix2(within_bounds, _half(bounds)), basis.monomials),
        )
    end
end

function half_newton_polytope(
    p::SA.AlgebraElement{<:MB.Algebra{BT,B,M}},
    gs::AbstractVector{<:SA.AlgebraElement},
    vars,
    maxdegree,
    newton::NewtonDegreeBounds,
) where {BT,B,M}
    if !is_commutative(vars)
        throw(
            ArgumentError(
                "Multipartite Newton polytope not supported with noncommutative variables.",
            ),
        )
    end
    parts = newton.variable_groups
    if !all(
        i -> all(j -> i == j || isempty(parts[i] ∩ parts[j]), eachindex(parts)),
        eachindex(parts),
    )
        throw(
            ArgumentError(
                "Parts are not disjoint in multipartite Newton polytope estimation: $parts.",
            ),
        )
    end
    # Some variables might be in no part...
    missing_vars = setdiff(vars, reduce(union, parts))
    if isempty(missing_vars)
        all_parts = parts
    else
        # in that case, create a part with the missing ones
        all_parts = (parts..., missing_vars)
    end
    if length(all_parts) == 1
        # all variables on same part, fallback to shortcut
        return half_newton_polytope(
            p,
            gs,
            vars,
            maxdegree,
            NewtonDegreeBounds(tuple()),
        )
    end
    bases = map(
        part -> half_newton_polytope(
            p,
            gs,
            part,
            maxdegree,
            NewtonDegreeBounds(tuple()),
        ),
        all_parts,
    )
    bounds = putinar_degree_bounds(p, gs, vars, maxdegree)
    return _cartesian_product([b[1] for b in bases], bounds),
    [
        _cartesian_product(
            [b[2][i] for b in bases],
            minus_shift(bounds, gs[i]),
        ) for i in eachindex(first(bases)[2])
    ]
end

function half_newton_polytope(
    p::SA.AlgebraElement{<:MB.Algebra{BT,B,M}},
    gs::AbstractVector{<:SA.AlgebraElement},
    vars,
    maxdegree,
    ::NewtonDegreeBounds{Tuple{}},
) where {BT,B,M}
    if is_commutative(vars)
        # TODO take `variable_groups` into account
        bounds = putinar_degree_bounds(p, gs, vars, maxdegree)
        full = MB.FullBasis{B,M}()
        return maxdegree_gram_basis(full, _half(bounds)),
        MB.explicit_basis_type(typeof(full))[
            multiplier_basis(g, bounds) for g in gs
        ]
    else
        if !isempty(gs)
            error(
                "Inequalities constraints not supported with noncommutative variables",
            )
        end
        # Non-commutative variables
        # We use Newton chip method of [Section 2.3, BKP16].
        #
        # [BKP16] Sabine Burgdorf, Igor Klep, and Janez Povh.
        # *Optimization of polynomials in non-commuting variables*.
        # Berlin: Springer, 2016.
        vars = unique!(sort(vars))
        bounds = _half(putinar_degree_bounds(p, gs, vars, maxdegree))
        monos = MP.monomial_type(typeof(p))[]
        for mono in MB.explicit_basis(p).monomials
            if _is_hermitian_square(mono)
                for i in 1:div(MP.degree(mono), 2)
                    w = _chip(mono, -i)
                    if within_bounds(w, bounds)
                        push!(monos, w)
                    end
                end
            end
        end
        basis = MB.SubBasis{B}(monos)
        return basis, typeof(basis)[]
    end
end

function half_newton_polytope(
    p::SA.AlgebraElement,
    gs::AbstractVector{<:MP.AbstractPolynomialLike},
    vars,
    maxdegree,
    filter,
)
    return half_newton_polytope(
        p,
        MB.algebra_element.(gs),
        vars,
        maxdegree,
        filter,
    )
end

function half_newton_polytope(
    p::SA.AlgebraElement,
    gs::AbstractVector{<:SA.AlgebraElement{<:MB.Algebra{B},T}},
    vars,
    maxdegree,
    filter::NewtonFilter{<:NewtonDegreeBounds},
) where {B,T}
    basis, multipliers_bases =
        half_newton_polytope(p, gs, vars, maxdegree, filter.outer_approximation)
    bases = copy(multipliers_bases)
    push!(bases, basis)
    gs = copy(gs)
    push!(gs, MB.constant_algebra_element(B, T))
    filtered_bases = post_filter(p, gs, bases)
    # The last one will be recomputed by the ideal certificate
    return filtered_bases[end], filtered_bases[1:(end-1)]
end

struct SignChange{T}
    sign::T
    Δ::Int
end
Base.copy(s::SignChange) = s
Base.iszero(::SignChange) = false
MA.scaling_convert(::Type, s::SignChange) = s
Base.:*(s::SignChange, α::Real) = SignChange(s.sign * α, s.Δ)

struct SignCount
    unknown::Int
    positive::Int
    negative::Int
end
SignCount() = SignCount(0, 0, 0)
Base.iszero(::SignCount) = false
function _sign(c::SignCount)
    if !iszero(c.unknown)
        return missing
    elseif iszero(c.positive)
        return -1
    elseif iszero(c.negative)
        return -1
    else
        return missing
    end
end

function Base.:*(α, a::SignCount)
    if α > 0
        return a
    elseif α < 0
        return SignCount(a.unknown, a.negative, a.positive)
    else
        error("Cannot multiply `SignCount`` with `$α`")
    end
end

Base.:*(a::SignCount, α) = α * a

function Base.:+(a::SignCount, b::SignCount)
    return SignCount(
        a.unknown + b.unknown,
        a.positive + b.positive,
        a.negative + b.negative,
    )
end

function Base.:+(c::SignCount, a::SignChange{Missing})
    #@assert c.unknown >= -a.Δ
    return SignCount(c.unknown + a.Δ, c.positive, c.negative)
end

function Base.:+(c::SignCount, a::SignChange{<:Number})
    if a.sign > 0
        #@assert c.positive >= -a.Δ
        return SignCount(c.unknown, c.positive + a.Δ, c.negative)
    elseif a.sign < 0
        #@assert c.negative >= -a.Δ
        return SignCount(c.unknown, c.positive, c.negative + a.Δ)
    elseif iszero(a.sign)
        error(
            "A polynomial should never contain a term with zero coefficient but found `$(a.sign)`.",
        )
    else
        error("Cannot determine sign of `$(a.sign)`.")
    end
end

Base.convert(::Type{SignCount}, Δ::SignChange) = SignCount() + Δ

function increase(cache, counter, generator_sign, monos, mult)
    for a in monos
        for b in monos
            MA.operate_to!(
                cache,
                *,
                MB.algebra_element(mult),
                MB.algebra_element(a),
                MB.algebra_element(b),
            )
            MA.operate!(
                SA.UnsafeAddMul(*),
                counter,
                _term_constant_monomial(
                    SignChange((a != b) ? missing : generator_sign, 1),
                    mult,
                ),
                cache,
            )
        end
    end
end

struct _DictCoefficients{K,V} <: SA.AbstractCoefficients{K,V}
    inner::Dict{K,V}
end

SA.nonzero_pairs(d::_DictCoefficients) = d.inner
Base.keys(d::_DictCoefficients) = keys(d.inner)
Base.haskey(d::_DictCoefficients, k) = haskey(d.inner, k)

Base.getindex(d::_DictCoefficients{K}, key::K) where {K} = d.inner[key]

function SA.unsafe_push!(c::_DictCoefficients{K}, key::K, value) where {K}
    c.inner[key] = if haskey(c.inner, key)
        MA.operate!!(+, c.inner[key], value)
    else
        value
    end
    return c
end

function _term(α, p::MB.Polynomial{B,M}) where {B,M}
    return MB.algebra_element(MP.term(α, p.monomial), MB.FullBasis{B,M}())
end

function _term_constant_monomial(α, ::MB.Polynomial{B,M}) where {B,M}
    return _term(α, MB.Polynomial{B}(MP.constant_monomial(M)))
end

# If `mono` is such that there is no other way to have `mono^2` by multiplying
# two different monomials of `monos` and `mono` is not in `X` then, the corresponding
# diagonal entry of the Gram matrix will be zero hence the whole column and row
# will be zero hence we can remove this monomial.
# See [Proposition 3.7, CLR95], [Theorem 2, L09] or [Section 2.4, BKP16].

# This is generalized here to the case with constraints as detailed in [L23].

# [CLR95] Choi, M. D. and Lam, T. Y. and Reznick, B.
# *Sum of Squares of Real Polynomials*.
# Proceedings of Symposia in Pure mathematics (1995)
#
# [L09] Löfberg, Johan.
# *Pre-and post-processing sum-of-squares programs in practice*.
# IEEE transactions on automatic control 54.5 (2009): 1007-1011.
#
# [BKP16] Sabine Burgdorf, Igor Klep, and Janez Povh.
# *Optimization of polynomials in non-commuting variables*.
# Berlin: Springer, 2016.
#
# [L23] Legat, Benoît
# *Exploiting the Structure of a Polynomial Optimization Problem*
# SIAM Conference on Applications of Dynamical Systems, 2023
function post_filter(
    poly::SA.AlgebraElement,
    generators,
    multipliers_gram_monos,
)
    # We use `_DictCoefficients` instead `SA.SparseCoefficients` because
    # we need to keep it canonicalized (without duplicate actually)
    # and don't care about the list of monomials being ordered
    counter = MB.algebra_element(
        _DictCoefficients(Dict{MP.monomial_type(typeof(poly)),SignCount}()),
        MB.implicit_basis(SA.basis(poly)),
    )
    algebra = MB.algebra(MB.implicit_basis(SA.basis(poly)))
    cache = zero(Float64, algebra)
    cache2 = zero(Float64, algebra)
    cache3 = zero(SignCount, algebra)
    cache4 = zero(SignCount, algebra)
    for (mono, v) in SA.nonzero_pairs(SA.coeffs(poly))
        MA.operate!(
            SA.UnsafeAddMul(*),
            counter,
            _term(SignChange(_sign(v), 1), SA.basis(poly)[mono]),
        )
    end
    for (mult, gram_monos) in zip(generators, multipliers_gram_monos)
        for (mono, v) in SA.nonzero_pairs(SA.coeffs(mult))
            increase(
                cache,
                counter,
                -_sign(v),
                gram_monos,
                SA.basis(mult)[mono],
            )
        end
    end
    function decrease(sign, a, b, generator)
        MA.operate_to!(
            cache,
            *,
            MB.algebra_element(a),
            MB.algebra_element(b),
            generator,
        )
        MA.operate_to!(
            cache3,
            *,
            _term(SignChange(sign, -1), a),
            MB.algebra_element(b),
        )
        MA.operate_to!(
            cache4,
            *,
            cache3,
            generator,
        )
        MA.operate!(
            SA.UnsafeAddMul(*),
            counter,
            cache4,
        )
        for mono in SA.supp(cache4)
            if !haskey(SA.coeffs(counter), SA.basis(counter)[mono])
                @show counter
                @show SA.basis(counter)[mono]
            end
            count = SA.coeffs(counter)[SA.basis(counter)[mono]]
            count_sign = _sign(count)
            # This means the `counter` has a sign and it didn't have a sign before
            # so we need to delete back edges
            if !ismissing(count_sign) &&
               (ismissing(count) || count != count_sign)
                # TODO could see later if deleting the counter improves perf
                if haskey(back, mono)
                    for (i, j) in back[mono]
                        delete(i, j)
                    end
                end
            end
        end
    end
    back = Dict{eltype(eltype(multipliers_gram_monos)),Vector{Tuple{Int,Int}}}()
    keep = [ones(Bool, length(monos)) for monos in multipliers_gram_monos]
    function delete(i, j)
        if !keep[i][j]
            return
        end
        keep[i][j] = false
        a = multipliers_gram_monos[i][j]
        decrease(-1, a, a, generators[i])
        for (k, b) in enumerate(multipliers_gram_monos[i])
            if keep[i][k]
                decrease(missing, a, b, generators[i])
                decrease(missing, b, a, generators[i])
            end
        end
    end
    for i in eachindex(generators)
        for (j, mono) in enumerate(multipliers_gram_monos[i])
            MA.operate_to!(
                cache,
                *,
                # Dummy coef to help convert to `SignCount` which is the `eltype` of `cache`
                _term(1.0, mono),
                MB.algebra_element(mono),
            )
            # The `eltype` of `cache` is `SignCount`
            # so there is no risk of term cancellation
            MA.operate_to!(
                cache2,
                *,
                cache,
                generators[i],
            )
            for w in SA.supp(cache2)
                if ismissing(
                    _sign(SA.coeffs(counter)[SA.basis(counter)[w]]),
                )
                    push!(get!(back, w, Tuple{Int,Int}[]), (i, j))
                else
                    delete(i, j)
                end
            end
        end
    end
    return [
        _sub(gram_monos, keep) for
        (keep, gram_monos) in zip(keep, multipliers_gram_monos)
    ]
end

function _sub(basis::SubBasis{B}, I) where {B}
    return SubBasis{B}(basis.monomials[I])
end

function _weight_type(::Type{T}, ::Type{BT}) where {T,BT}
    return SA.AlgebraElement{
        MA.promote_operation(
            MB.algebra,
            MA.promote_operation(MB.implicit_basis, BT),
        ),
        T,
        MA.promote_operation(
            MB.sparse_coefficients,
            MP.polynomial_type(MP.monomial_type(BT), T),
        ),
    }
end

function half_newton_polytope(a::SA.AlgebraElement, vars, filter)
    return half_newton_polytope(
        a,
        _weight_type(Bool, typeof(SA.basis(a)))[],
        vars,
        _maxdegree(a, vars),
        filter,
    )[1]
end

function half_newton_polytope(a::SA.AlgebraElement, filter)
    return half_newton_polytope(a, MP.variables(a), filter)
end

function half_newton_polytope(basis::MB.SubBasis, args...)
    a = MB.algebra_element(
        SA.SparseCoefficients(basis.monomials, ones(length(basis))),
        MB.implicit_basis(basis),
    )
    return half_newton_polytope(a, args...)
end

function monomials_half_newton_polytope(monos, args...)
    return half_newton_polytope(
        MB.SubBasis{MB.Monomial}(monos),
        args...,
    ).monomials
end
