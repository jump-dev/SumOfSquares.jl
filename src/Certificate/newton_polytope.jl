cfld(x::NTuple{2,Int}, n) = (cld(x[1], n), fld(x[2], n))

function sub_extdegree(X::AbstractVector{<:MP.AbstractMonomial}, vars)
    if isempty(X)
        return (0, 0)
    else
        return extrema(map(mono -> sum(var -> MP.degree(mono, var), vars), X))
    end
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

function _filter(X::AbstractVector{<:MP.AbstractMonomial}, extdeg, exp, n)
    mindeg, maxdeg = cfld(extdeg, 2)
    minmultideg, maxmultideg = Vector{Int}(undef, n), Vector{Int}(undef, n)
    for i in 1:n
        exponent_i(mono) = exp(mono, i)
        minmultideg[i] = cld(mapreduce(exponent_i, min, X), 2)
        maxmultideg[i] = fld(mapreduce(exponent_i, max, X), 2)
    end
    return mindeg,
    maxdeg,
    mono -> begin
        all(i -> minmultideg[i] <= exp(mono, i) <= maxmultideg[i], 1:n)
    end
end
function _full_filter(X::AbstractVector{<:MP.AbstractMonomial}, extdeg, exp, n)
    mindeg, maxdeg, filter = _filter(X, extdeg, exp, n)
    return mono -> mindeg <= MP.degree(mono) <= maxdeg && filter(mono)
end

function _sub_half_newton_polytope(
    X::AbstractVector{<:MP.AbstractMonomial},
    extdeg,
    exp,
    vars,
)
    mindeg, maxdeg, filter = _filter(X, extdeg, exp, length(vars))
    return MP.monomials(vars, mindeg:maxdeg, filter)
end

function sub_half_newton_polytope(
    X::AbstractVector{<:MP.AbstractMonomial},
    vars,
)
    return _sub_half_newton_polytope(
        X,
        sub_extdegree(X, vars),
        (mono, i) -> MP.degree(mono, vars[i]),
        vars,
    )
end

function is_commutative(vars)
    return length(vars) < 2 || prod(vars[1:2]) == prod(reverse(vars[1:2]))
end

abstract type AbstractNewtonPolytopeApproximation end

# `maxdegree` and `mindegree` on each variable separately
# and on each group of variable.
struct NewtonDegreeBounds{NPT} <: AbstractNewtonPolytopeApproximation
    variable_groups::NPT
end

# Multipartite
# TODO we might do this recursively : do 2 parts, merge them, merge with next
#      one and so on so that the filter at the end prunes more.
function half_newton_polytope(X::AbstractVector, newton::NewtonDegreeBounds)
    if !is_commutative(MP.variables(X))
        throw(
            ArgumentError(
                "Multipartite Newton polytope not supported with noncommutative variables.",
            ),
        )
    end
    parts = newton.variable_groups
    if !all(
        i -> all(j -> i == j || isempty(parts[i] ∩ parts[j]), 1:length(parts)),
        1:length(parts),
    )
        throw(
            ArgumentError(
                "Parts are not disjoint in multipartite Newton polytope estimation: $parts.",
            ),
        )
    end
    # Some variables might be in no part...
    missing = setdiff(MP.variables(X), reduce(union, parts))
    if isempty(missing)
        all_parts = parts
    else
        # in that case, create a part with the missing ones
        all_parts = (parts..., missing)
    end
    if length(all_parts) == 1
        # all variables on same part, fallback to shortcut
        return half_newton_polytope(X, NewtonDegreeBounds(tuple()))
    end
    monomial_vectors = map(vars -> sub_half_newton_polytope(X, vars), all_parts)
    # Cartesian product of the newton polytopes of the different parts
    product = [prod(monos) for monos in Iterators.product(monomial_vectors...)]
    mindeg, maxdeg = cfld(MP.extdegree(X), 2)
    # We know that the degree inequalities are satisfied variable-wise and
    # part-wise but for all variables together so we filter with that
    gram_monos = filter(mono -> mindeg <= MP.degree(mono) <= maxdeg, product)
    return MP.monomial_vector(gram_monos)
end

# Filters out points ouside the Newton polytope from the
# outer approximation given by `outer_approximation`.
struct NewtonFilter{N<:AbstractNewtonPolytopeApproximation} <:
       AbstractNewtonPolytopeApproximation
    outer_approximation::N
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

# Shortcut for more efficient `extdeg` and `exp` function in case all the
# variables are in the same part
function half_newton_polytope(X::AbstractVector, ::NewtonDegreeBounds{Tuple{}})
    vars = MP.variables(X)
    if is_commutative(vars)
        # Commutative variables
        exp(mono, i) = MP.exponents(mono)[i]
        return _sub_half_newton_polytope(X, MP.extdegree(X), exp, vars)
    else
        # Non-commutative variables
        # We use Newton chip method of [Section 2.3, BKP16].
        #
        # [BKP16] Sabine Burgdorf, Igor Klep, and Janez Povh.
        # *Optimization of polynomials in non-commuting variables*.
        # Berlin: Springer, 2016.
        vars = unique!(sort(vars))
        function ncexp(mono, i)
            mvars = MP.variables(mono)
            return mapreduce(
                j -> mvars[j] == vars[i] ? MP.exponents(mono)[j] : 0,
                +,
                eachindex(mvars),
                init = 0,
            )
        end
        filter = _full_filter(X, MP.extdegree(X), ncexp, length(vars))
        _monos = eltype(X)[]
        for mono in X
            if _is_hermitian_square(mono)
                for i in 1:div(MP.degree(mono), 2)
                    w = _chip(mono, -i)
                    if filter(w)
                        push!(_monos, w)
                    end
                end
            end
        end
        return MP.monomial_vector(_monos)
    end
end

function half_newton_polytope(basis::SA.ExplicitBasis, newton::NewtonFilter)
    gram_monos = half_newton_polytope(monos, newton.outer_approximation)
    return post_filter(gram_monos, monos)
end

# If `mono` is such that there is no other way to have `mono^2` by multiplying
# two different monomials of `monos` and `mono` is not in `X` then, the corresponding
# diagonal entry of the Gram matrix will be zero hence the whole column and row
# will be zero hence we can remove this monomial.
# See [Proposition 3.7, CLR95], [Theorem 2, L09] or [Section 2.4, BKP16].

# [CLR95] Choi, M. D. and Lam, T. Y. and Reznick, B.
# *Sum of Squares of Real Polynomials*.
# Proceedings of Symposia in Pure mathematics (1995)
#
# [L09] Lofberg, Johan.
# *Pre-and post-processing sum-of-squares programs in practice*.
# IEEE transactions on automatic control 54.5 (2009): 1007-1011.
#
# [BKP16] Sabine Burgdorf, Igor Klep, and Janez Povh.
# *Optimization of polynomials in non-commuting variables*.
# Berlin: Springer, 2016.
function post_filter(monos, X)
    num = Dict(mono => 1 for mono in X)
    function _increase(mono)
        return num[mono] = get(num, mono, 0) + 1
    end
    function _decrease(mono)
        value = num[mono] - 1
        if iszero(value)
            delete!(num, mono)
        else
            num[mono] = value
        end
        return iszero(value)
    end
    for a in monos
        for b in monos
            if a != b
                _increase(a * b)
            end
        end
    end
    back = Dict{eltype(monos),Int}()
    keep = ones(Bool, length(monos))
    function _delete(i)
        keep[i] = false
        a = monos[i]
        for (j, b) in enumerate(monos)
            if keep[j] && a != b
                for w in [a * b, b * a]
                    if _decrease(w) && haskey(back, w)
                        _delete(back[w])
                    end
                end
            end
        end
    end
    for (i, mono) in enumerate(monos)
        w = mono^2
        if haskey(num, w)
            back[w] = i
        else
            _delete(i)
        end
    end
    return monos[findall(keep)]
end

struct DegreeBounds{M}
    mindegree::Int
    maxdegree::Int
    variablewise_mindegree::M
    variablewise_maxdegree::M
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
min_degree(p::MB.Polynomial{MB.Monomial}, v) = MP.degree(p.monomial, v)
minus_min_degree(p, v) = -min_degree(p, v)

function _monomial(vars, exps)
    if any(Base.Fix2(isless, 0), exps)
        return
    end
    return prod(vars .^ exps)
end

function max_degree(p::SA.AlgebraElement, v)
    return mapreduce(Base.Fix2(max_degree, v), max, SA.supp(p))
end
max_degree(p::MB.Polynomial{MB.Monomial}, v) = MP.degree(p.monomial, v)

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
    return DegreeBounds(
        min_shift(d.mindegree, _mindegree(p)),
        d.maxdegree - _maxdegree(p),
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
_sign(a) = missing

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

function deg_range(deg, p, gs, gram_deg, truncation)
    d_max = min(deg(p), truncation)
    sign = deg_sign(deg, p, d_max)
    for g in gs
        d_g, sign_g = deg_sign(deg, g)
        d_s = gram_deg(g)
        if isempty(d_s) || d_g + 2minimum(d_s) > truncation
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
# Maximum value of `deg(s_0 = p - sum s_i g_i for g in gs) in range` where
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
_mindegree(a::SA.AlgebraElement) = minimum(_mindegree, SA.supp(a))
_maxdegree(a::SA.AlgebraElement) = maximum(_maxdegree, SA.supp(a))
#_mindegree(basis::MB.SubBasis{MB.Monomial}) = MP.mindegree(basis.monomials)
#_maxdegree(basis::MB.SubBasis{MB.Monomial}) = MP.maxdegree(basis.monomials)
_mindegree(p::MB.Polynomial{MB.Monomial}) = MP.mindegree(p.monomial)
_maxdegree(p::MB.Polynomial{MB.Monomial}) = MP.maxdegree(p.monomial)

function putinar_degree_bounds(
    p::SA.AlgebraElement,
    gs::AbstractVector{<:SA.AlgebraElement},
    vars,
    maxdegree,
)
    mindegree = 0
    # TODO homogeneous case
    mindeg(g) = _min_half(min_shift(mindegree, _mindegree(g)))
    maxdeg(g) = _max_half(maxdegree - _maxdegree(g))
    degrange(g) = mindeg(g):maxdeg(g)
    minus_degrange(g) = -maxdeg(g):-mindeg(g)
    # The multiplier will have degree `0:2fld(maxdegree - _maxdegree(g), 2)`
    mindegree =
        -deg_range(p -> -_mindegree(p), p, gs, minus_degrange, -maxdegree:0)
    if isnothing(mindegree)
        return
    end
    maxdegree = deg_range(_maxdegree, p, gs, degrange, 0:maxdegree)
    if isnothing(maxdegree)
        return
    end
    vars_mindeg = map(vars) do v
        return -deg_range(
            Base.Fix2(minus_min_degree, v),
            p,
            gs,
            minus_degrange,
            -maxdegree:0,
        )
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

function multiplier_basis(g::SA.AlgebraElement, bounds::DegreeBounds)
    shifted = minus_shift(bounds, g)
    if isnothing(shifted)
        halved = nothing
    else
        halved = _half(shifted)
    end
    basis = MB.FullBasis{MB.Monomial,MP.monomial_type(typeof(g))}()
    if isnothing(halved)
        # TODO add `MB.empty_basis` to API
        return MB.maxdegree_basis(
            basis,
            MP.variables(bounds.variablewise_mindegree),
            -1,
        )
    else
        return maxdegree_gram_basis(basis, halved)
    end
end

function half_newton_polytope(
    p::MP.AbstractPolynomialLike,
    gs::AbstractVector{<:MP.AbstractPolynomialLike},
    vars,
    maxdegree,
    ::NewtonDegreeBounds,
)
    # TODO take `variable_groups` into account
    bounds = putinar_degree_bounds(p, gs, vars, maxdegree)
    return [multiplier_basis(g, bounds) for g in gs]
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
        [MB.algebra_element(MP.coefficients(g), MB.SubBasis{MB.Monomial}(MP.monomials(g))) for g in gs],
        vars,
        maxdegree,
        filter,
    )
end

function half_newton_polytope(
    p::SA.AlgebraElement,
    gs::AbstractVector{<:SA.AlgebraElement},
    vars,
    maxdegree,
    ::NewtonFilter{<:NewtonDegreeBounds},
)
    bounds = putinar_degree_bounds(p, gs, vars, maxdegree)
    bases = [multiplier_basis(g, bounds) for g in gs]
    push!(
        bases,
        maxdegree_gram_basis(
            MB.FullBasis{MB.Monomial,MP.monomial_type(typeof(p))}(),
            _half(bounds),
        ),
    )
    gs = copy(gs)
    push!(gs, MB.constant_algebra_element(MA.promote_operation(SA.basis, eltype(gs)), eltype(eltype(gs))))
    filtered_bases = post_filter(p, gs, bases)
    # The last one will be recomputed by the ideal certificate
    return filtered_bases[1:(end-1)]
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

function Base.:+(c::SignCount, a::SignChange{Missing})
    @assert c.unknown >= -a.Δ
    return SignCount(c.unknown + a.Δ, c.positive, c.negative)
end
function Base.:+(c::SignCount, a::SignChange{<:Number})
    if a.sign > 0
        @assert c.positive >= -a.Δ
        return SignCount(c.unknown, c.positive + a.Δ, c.negative)
    elseif a.sign < 0
        @assert c.negative >= -a.Δ
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
                SignChange((a != b) ? missing : generator_sign, 1,),
                cache,
            )
        end
    end
end

function post_filter(poly::SA.AlgebraElement, generators, multipliers_gram_monos)
    counter = MB.algebra_element(
        zero(MP.polynomial_type(MP.monomial_type(typeof(SA.basis(poly))), SignCount)),
        MB.implicit_basis(SA.basis(poly)),
    )
    cache = MB.algebra_element(
        zero(MP.polynomial_type(MP.monomial_type(typeof(SA.basis(poly))), Float64)),
        MB.implicit_basis(SA.basis(poly)),
    )
    for (mono, v) in SA.nonzero_pairs(poly)
        MA.operate!(
            SA.UnsafeAddMul(*),
            counter,
            SignChange(_sign(v), 1),
            mono,
        )
    end
    for (mult, gram_monos) in zip(generators, multipliers_gram_monos)
        for (mono, v) in SA.nonzero_pairs(mult)
            increase(cache, counter, -_sign(v), gram_monos, mono)
        end
    end
    function decrease(sign, a, b, c)
        MA.operate_to!(
            cache,
            *,
            MB.algebra_element(a),
            MB.algebra_element(b),
            MB.algebra_element(c),
        )
        MA.operate!(SA.UnsafeAddMul(*), counter, SignChange(sign, -1), cache)
        MA.operate!(SA.canonical, counter)
        for mono in SA.supp(cache)
            count = counter[mono]
            count_sign = _sign(count)
            # This means the `counter` has a sign and it didn't have a sign before
            # so we need to delete back edges
            if !ismissing(count_sign) && (ismissing(count) || count != count_sign)
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
        for (k, v) in SA.nonzero_pairs(generators[i])
            sign = -_sign(v)
            decrease(sign, k, a, a)
            for (j, b) in enumerate(multipliers_gram_monos[i])
                if keep[i][j]
                    decrease(missing, k, a, b)
                    decrease(missing, k, b, a)
                end
            end
        end
    end
    for i in eachindex(generators)
        for k in SA.supp(generators[i])
            for (j, mono) in enumerate(multipliers_gram_monos[i])
                MA.operate_to!(
                    cache,
                    *,
                    MB.algebra_element(k),
                    MB.algebra_element(mono),
                    MB.algebra_element(mono),
                )
                for w in SA.supp(cache)
                    if ismissing(_sign(SA.coeffs(counter)[SA.basis(counter)[w]]))
                        push!(get(back, w, Tuple{Int,Int}[]), (i, j))
                    else
                        delete(i, j)
                    end
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

function half_newton_polytope(a::SA.AlgebraElement, args...)
    @assert SA.basis(a) isa MB.SubBasis{MB.Monomial}
    half_newton_polytope(SA.coeffs(a, MB.FullBasis{MB.Monomial,MP.monomial_type(typeof(a))}()), args...)
end
