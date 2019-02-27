# Inspired from SOSTools
export randpsd, randsos

cfld(x::NTuple{2,Int}, n) = (cld(x[1], n), fld(x[2], n))

function sub_extdegree(X::AbstractVector{<:MP.AbstractMonomial}, vars)
    if isempty(X)
        return (0, 0)
    else
        return extrema(map(mono -> sum(var -> degree(mono, var), vars), X))
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

function _sub_half_newton_polytope(X::AbstractVector{<:MP.AbstractMonomial},
                                   extdeg, exp, vars)
    mindeg, maxdeg = cfld(extdeg, 2)
    n = length(vars)
    minmultideg, maxmultideg = Vector{Int}(undef, n), Vector{Int}(undef, n)
    for i in 1:n
        exponent_i(mono) = exp(mono, i)
        minmultideg[i] = cld(mapreduce(exponent_i, min, X), 2)
        maxmultideg[i] = fld(mapreduce(exponent_i, max, X), 2)
    end
    monomials(vars, mindeg:maxdeg,
              mono -> all(i -> minmultideg[i] <= exp(mono, i) <= maxmultideg[i],
                          1:n))
end

function sub_half_newton_polytope(X::AbstractVector{<:MP.AbstractMonomial},
                                  vars)
    _sub_half_newton_polytope(X, sub_extdegree(X, vars),
                              (mono, i) -> degree(mono, vars[i]), vars)
end

# Multipartite
# TODO we might do this recursively : do 2 parts, merge them, merge with next
#      one and so on so that the filter at the end prunes more.
function half_newton_polytope(X::AbstractVector, parts::Tuple)
    if !all(i -> all(j -> i == j || isempty(parts[i] ∩ parts[j]),
                     1:length(parts)),
            1:length(parts))
        throw(ArgumentError("Parts are not disjoint in multipartite Newton polytope estimation: $parts"))
    end
    # Some variables might be in no part...
    missing = setdiff(variables(X), reduce(union, parts))
    if isempty(missing)
        all_parts = parts
    else
        # in that case, create a part with the missing ones
        all_parts = (parts..., missing)
    end
    if length(all_parts) == 1
        # all variables on same part, fallback to shortcut
        return half_newton_polytope(X, tuple())
    end
    monovecs = map(vars -> sub_half_newton_polytope(X, vars), all_parts)
    # Cartesian product of the newton polytopes of the different parts
    product = [prod(monos) for monos in Iterators.product(monovecs...)]
    mindeg, maxdeg = cfld(extdegree(X), 2)
    # We know that the degree inequalities are satisfied variable-wise and
    # part-wise but for all variables together so we filter with that
    return monovec(filter(mono -> mindeg <= degree(mono) <= maxdeg, product))
end
# Shortcut for more efficient `extdeg` and `exp` function in case all the
# variables are in the same part
function half_newton_polytope(X::AbstractVector, parts::Tuple{})
    return _sub_half_newton_polytope(X, extdegree(X),
                                     (mono, i) -> exponents(mono)[i],
                                     variables(X))
end

function monomials_half_newton_polytope(X::AbstractVector, parts)
    half_newton_polytope(monovec(X), parts)
end

function randpsd(n; r=n, eps=0.1)
    Q = randn(n,n)
    d = zeros(Float64, n)
    d[1:r] = eps .+ abs.(randn(r))
    Q' * Diagonal(d) * Q
end

function _randsos(X::AbstractVector{<:MP.AbstractMonomial}; r=-1, monotype=:Classic, eps=0.1)
    if monotype == :Classic
        x = monomials_half_newton_polytope(X, tuple())
    elseif monotype == :Gram
        x = X
    else
        throw(ArgumentError("Monotype $monotype not known"))
    end
    n = length(x)
    if r < 0
        r = n
    end
    GramMatrix(randpsd(n, r=r, eps=eps), x)
end

randsos(X::AbstractVector; kws...) = _randsos(monovec(X); kws...)
