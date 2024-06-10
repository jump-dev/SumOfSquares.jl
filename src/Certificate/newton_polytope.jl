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
_mindegree(a::SA.AlgebraElement) = minimum(_mindegree, SA.supp(a))
_maxdegree(a::SA.AlgebraElement) = maximum(_maxdegree, SA.supp(a), init = 0)
#_mindegree(basis::MB.SubBasis{MB.Monomial}) = MP.mindegree(basis.monomials)
#_maxdegree(basis::MB.SubBasis{MB.Monomial}) = MP.maxdegree(basis.monomials)
function _mindegree(p::MB.Polynomial{B}) where {B}
    if _is_monomial_basis(B)
        MP.mindegree(p.monomial)
    else
        error("TODO $B")
    end
end
_maxdegree(p::MB.Polynomial) = MP.maxdegree(p.monomial)

_is_monomial_basis(::Type{<:MB.AbstractMonomialIndexed}) = false
_is_monomial_basis(::Type{<:Union{MB.Monomial,MB.ScaledMonomial}}) = true
_is_monomial_basis(::Type{<:SA.AlgebraElement{<:MB.Algebra{BT,B}}}) where {BT,B} = _is_monomial_basis(B)

# Minimum degree of a gram basis for a gram matrix `s`
# such that the minimum degree of `s * g` is at least `mindegree`.
function _multiplier_mindegree(
    mindegree,
    g::SA.AlgebraElement,
)
    if _is_monomial_basis(typeof(g))
        return _min_half(min_shift(mindegree, _mindegree(g)))
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
function _multiplier_maxdegree(maxdegree, g::SA.AlgebraElement)
    return _max_half(maxdegree - _maxdegree(g))
end

function _multiplier_deg_range(range, g::SA.AlgebraElement)
    return _multiplier_mindegree(minimum(range), g):_multiplier_maxdegree(maximum(range), g)
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
    degrange(g) = _multiplier_deg_range(mindegree:maxdegree, g)
    minus_degrange(g) = (-).(degrange(g))
    # The multiplier will have degree `0:2fld(maxdegree - _maxdegree(g), 2)`
    mindegree = if _is_monomial_basis(typeof(p))
        -deg_range((-) ∘ _mindegree, p, gs, minus_degrange, -maxdegree:0)
    else
        0
    end
    if isnothing(mindegree)
        return
    end
    maxdegree = deg_range(_maxdegree, p, gs, degrange, 0:maxdegree)
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

function multiplier_basis(g::SA.AlgebraElement{<:MB.Algebra{BT,B}}, bounds::DegreeBounds) where {BT,B}
    shifted = minus_shift(bounds, g)
    if isnothing(shifted)
        halved = nothing
    else
        halved = _half(shifted)
    end
    basis = MB.FullBasis{B,MP.monomial_type(typeof(g))}()
    if isnothing(halved)
        return MB.empty_basis(MB.explicit_basis_type(typeof(basis)))
    else
        return maxdegree_gram_basis(basis, halved)
    end
end

function half_newton_polytope(
    p::SA.AlgebraElement{<:MB.Algebra{BT,B,M}},
    gs::AbstractVector{<:SA.AlgebraElement},
    vars,
    maxdegree,
    ::NewtonDegreeBounds,
) where {BT,B,M}
    # TODO take `variable_groups` into account
    bounds = putinar_degree_bounds(p, gs, vars, maxdegree)
    full = MB.FullBasis{B,M}()
    return maxdegree_gram_basis(full, _half(bounds)),
        MB.explicit_basis_type(typeof(full))[multiplier_basis(g, bounds) for g in gs]
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
    filter::NewtonFilter{<:NewtonDegreeBounds},
)
    basis, multipliers_bases = half_newton_polytope(p, gs, vars, maxdegree, filter.outer_approximation)
    bases = copy(multipliers_bases)
    push!(bases, basis)
    gs = copy(gs)
    push!(gs, MB.constant_algebra_element(MA.promote_operation(SA.basis, eltype(gs)), eltype(eltype(gs))))
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

function Base.:+(a::SignCount, b::SignCount)
    return SignCount(a.unknown + b.unknown, a.positive + b.positive, a.negative + b.negative)
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

struct _DictCoefficients{K,V} <: SA.AbstractCoefficients{K,V}
    inner::Dict{K,V}
end

SA.nonzero_pairs(d::_DictCoefficients) = d.inner
Base.keys(d::_DictCoefficients) = keys(d.inner)

Base.getindex(d::_DictCoefficients{K}, key::K) where {K} = d.inner[key]

function SA.unsafe_push!(c::_DictCoefficients{K}, key::K, value) where {K}
    c.inner[key] = if haskey(c.inner, key)
        MA.operate!!(+, c.inner[key], value)
    else
        value
    end
    return c
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
function post_filter(poly::SA.AlgebraElement, generators, multipliers_gram_monos)
    # We use `_DictCoefficients` instead `SA.SparseCoefficients` because
    # we need to keep it canonicalized (without duplicate actually)
    # and don't care about the list of monomials being ordered
    counter = MB.algebra_element(_DictCoefficients(Dict{MP.monomial_type(typeof(poly)),SignCount}()), MB.implicit_basis(SA.basis(poly)))
    cache = zero(Float64, SA.algebra(MB.implicit_basis(SA.basis(poly))))
    for (mono, v) in SA.nonzero_pairs(poly)
        MA.operate!(
            SA.UnsafeAddMul(*),
            counter,
            SignChange(_sign(v), 1),
            MB.algebra_element(mono),
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
        for mono in SA.supp(cache)
            count = SA.coeffs(counter)[SA.basis(counter)[mono]]
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
        for (k, v) in SA.nonzero_pairs(SA.coeffs(generators[i]))
            mono = SA.basis(generators[i])[k]
            sign = -_sign(v)
            decrease(sign, mono, a, a)
            for (j, b) in enumerate(multipliers_gram_monos[i])
                if keep[i][j]
                    decrease(missing, mono, a, b)
                    decrease(missing, mono, b, a)
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
