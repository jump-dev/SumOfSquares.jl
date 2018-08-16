# Inspired from SOSTools
export getmonomialsforcertificate, randpsd, randsos

function map_extrema(f::Function, itr::AbstractVector)
    state = start(itr)
    it1, state = next(itr, state)
    mini = maxi = f(it1)
    while !done(itr, state)
        it, state = next(itr, state)
        fit = f(it)
        if fit < mini
            mini = fit
        end
        if fit > maxi
            maxi = fit
        end
    end
    mini, maxi
end

cfld(x::NTuple{2,Int}, n) = (cld(x[1], n), fld(x[2], n))

# TODO sparse with Newton polytope (Polyhedra.jl for convex hull)
function _getmonomialsforcertificate(X::AbstractVector{<:AbstractMonomial}, sparse=:No)
    if sparse == :No
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
        mindeg, maxdeg = cfld(extdegree(X), 2)
        n = nvariables(X)
        minmultideg, maxmultideg = Vector{Int}(undef, n), Vector{Int}(undef, n)
        for i in 1:n
            minmultideg[i], maxmultideg[i] = cfld(map_extrema(m -> exponents(m)[i], X), 2)
        end
        monomials(variables(X), mindeg:maxdeg,
                  m -> all(minmultideg .<= exponents(m) .<= maxmultideg))
    else
        error("Not supported yet :(")
    end
end
getmonomialsforcertificate(X::AbstractVector, sparse=:No) = _getmonomialsforcertificate(monovec(X), sparse)

function randpsd(n; r=n, eps=0.1)
    Q = randn(n,n)
    d = zeros(Float64, n)
    d[1:r] = eps .+ abs.(randn(r))
    Q' * Diagonal(d) * Q
end

function _randsos(X::AbstractVector{<:AbstractMonomial}; r=-1, monotype=:Classic, eps=0.1)
    if monotype == :Classic
        x = getmonomialsforcertificate(X)
    elseif monotype == :Gram
        x = X
    else
        throw(ArgumentError("Monotype $monotype not known"))
    end
    n = length(x)
    if r < 0
        r = n
    end
    MatPolynomial(randpsd(n, r=r, eps=eps), x)
end

randsos(X::AbstractVector; kws...) = _randsos(monovec(X); kws...)
