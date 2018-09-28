# Inspired from SOSTools
export getmonomialsforcertificate, randpsd, randsos

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
            exponent_i(m) = exponents(m)[i]
            minmultideg[i] = cld(Compat.mapreduce(exponent_i, min, X), 2)
            maxmultideg[i] = fld(Compat.mapreduce(exponent_i, max, X), 2)
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
    Q' * Compat.LinearAlgebra.Diagonal(d) * Q
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
