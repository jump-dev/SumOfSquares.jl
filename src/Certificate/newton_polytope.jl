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

# Multipartite
# TODO we might do this recursively : do 2 parts, merge them, merge with next
#      one and so on so that the filter at the end prunes more.
function half_newton_polytope(
    X::AbstractVector,
    parts::Tuple;
    apply_post_filter = true,
)
    if !is_commutative(MP.variables(X))
        throw(
            ArgumentError(
                "Multipartite Newton polytope not supported with noncommutative variables.",
            ),
        )
    end
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
        return half_newton_polytope(
            X,
            tuple();
            apply_post_filter = apply_post_filter,
        )
    end
    monovecs = map(vars -> sub_half_newton_polytope(X, vars), all_parts)
    # Cartesian product of the newton polytopes of the different parts
    product = [prod(monos) for monos in Iterators.product(monovecs...)]
    mindeg, maxdeg = cfld(MP.extdegree(X), 2)
    # We know that the degree inequalities are satisfied variable-wise and
    # part-wise but for all variables together so we filter with that
    monos =
        MP.monovec(filter(mono -> mindeg <= MP.degree(mono) <= maxdeg, product))
    if apply_post_filter
        return post_filter(monos, X)
    else
        return monos
    end
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
    return __chip(MP.constantmonomial(mono), 1, vars, exps, abs(n), op)
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
function half_newton_polytope(
    X::AbstractVector,
    ::Tuple{};
    apply_post_filter = true,
)
    vars = MP.variables(X)
    if is_commutative(vars)
        # Commutative variables
        exp(mono, i) = MP.exponents(mono)[i]
        monos = _sub_half_newton_polytope(X, MP.extdegree(X), exp, vars)
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
        monos = MP.monovec(_monos)
    end
    if apply_post_filter
        return post_filter(monos, X)
    else
        return monos
    end
end

# If `mono` is such that there is no other way to have `mono^2` by multiplying
# two different monomials of `monos` and `mono` is not in `X` then, the corresponding
# diagonal entry of the Gram matrix will be zero hence the whole column and row
# will be zero hence we can remove this monomial.
# See [Theorem 2, L09] or [Section 2.4, BKP16].
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

function monomials_half_newton_polytope(
    X::AbstractVector,
    parts;
    apply_post_filter = true,
)
    return half_newton_polytope(
        MP.monovec(X),
        parts;
        apply_post_filter = apply_post_filter,
    )
end

struct DegreeBounds{M}
    mindegree::Int
    maxdegree::Int
    variablewise_mindegree::M
    variablewise_maxdegree::M
end

function min_degree(p, v)
    return mapreduce(Base.Fix2(MP.degree, v), min, MP.monomials(p))
end
minus_min_degree(p, v) = -min_degree(p, v)

function _monomial(vars, exps)
    return prod(vars .^ exps)
end

function max_degree(p, v)
    return mapreduce(Base.Fix2(MP.degree, v), max, MP.monomials(p))
end

function minus_shift(
    deg,
    m::MP.AbstractMonomial,
    p::MP.AbstractPolynomialLike,
    offset,
)
    return _monomial(MP.variables(m), map(MP.powers(m)) do (v, d)
        return div(d - deg(p, v) + offset, 2)
    end)
end

function minus_shift(d::DegreeBounds, p::MP.AbstractPolynomialLike)
    return DegreeBounds(
        div(d - MP.mindegree(p) + 1, 2),
        div(d - MP.maxdegree(p), 2),
        minus_shift(min_degree, m, p, 1),
        minus_shift(min_degree, m, p, 0),
    )
end

_combine_sign(a, b) = (a == b ? a : zero(a))

function deg_range(deg, p, gs, gram_deg)
    d_max, sign = deg(p)
    for g in gs
        d_g, sign_g = deg(g) + 2gram_deg(g)
        # Multiply by `-1` because we move it to lhs
        # p = s_0 + sum s_i g_i -> p - sum s_i g_i = s_0
        sign_g = -sign_g
        if d_max == d_g
            sign = _combine_sign(sign, sign_g)
        end
        if d_max <= d_g
            d_max = d_g
        end
    end
    if iszero(sign) || (iseven(d_max) && sign == 1)
        return d_max
    else
        return d_max - 1
    end
end

function putinar_degree_bounds(
    p::MP.AbstractPolynialLike,
    gs::AbstractVector{<:MP.AbstractPolynialLike},
    vars,
    maxdegree,
)
    mindegree = 0
    # TODO homogeneous case
    minus_mindeg(g) = -max(0, div(mindegree - MP.mindegree(g) + 1, 2))
    # The multiplier will have degree `0:2div(maxdegree - MP.maxdegree(g), 2)`
    mindegree = -deg_range(p -> -MP.mindegree(p), p, gs, minus_mindeg)
    maxdeg(g) = div(maxdegree - MP.maxdegree(g), 2)
    maxdegree = deg_range(MP.maxdegree, p, gs, maxdeg)
    deg_range(MP.maxdegree, p, gs, maxdeg)
    vars_mindeg = map(vars) do v
        return -deg_range(Base.Fix2(minus_min_degree, v), p, gs, minus_mindeg)
    end
    vars_maxdeg = map(vars) do v
        return deg_range(Base.Fix2(max_degree, v), p, gs, maxdeg)
    end
    return DegreeBounds(
        mindegree,
        maxdegree,
        _monomial(vars, vars_mindeg),
        _monomial(vars, vars_maxdeg),
    )
end
