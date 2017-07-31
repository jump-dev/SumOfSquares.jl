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

# Cheap approximation of the convex hull as the approximation of:
#
# z such that mindeg < sum(z) < maxdeg
# |\
# |#\ <-------- sum(z) = maxdeg
# |##\
# |\<-\-------- sum(z) = mindeg
# | \##\
# +---------
cfld(x::NTuple{2,Int}, n) = (cld(x[1], 2), fld(x[2], 2))
# and:
#
# z such that minmultideg < z < maxmultideg
# | +----+ <--- maxmultideg
# | |####|
# | |####|
# | +----+
# | ^---------- minmultideg
# +---------
struct CheapOuterLibrary <: PolyhedraLibrary end
struct CheapOuterPolytope{N, T} <: Polyhedron{N, T}
    minmultideg::Vector{Int}
    maxmultideg::Vector{Int}
end
function newtonpolytope(Z::MonomialVector, lib::CheapOuterLibrary)
    n = length(Z.Z[1])
    minmultideg, maxmultideg = Vector{Int}(n), Vector{Int}(n)
    for i in 1:n
        a, b = map_extrema(z->z[i], Z.Z)
        minmultideg[i], maxmultideg[i] = cfld(map_extrema(z->z[i], Z.Z), 2)
    end
    CheapOuterPolytope{n, Int}(minmultideg, maxmultideg)
end
function Base.in(z::Vector{Int}, p::CheapOuterPolytope)
    reduce(&, true, p.minmultideg .<= z .<= p.maxmultideg)
end

function newtonpolytope(Z::MonomialVector, lib::PolyhedraLibrary)
    V = reduce(vcat, Matrix{Int}(0, 2), Z.Z.')
    vr = SimpleVRepresentation(V)
    polyhedron(vr, lib)
end

# TODO sparse with Newton polytope (Polyhedra.jl for convex hull)
function getmonomialsforcertificate(Z::MonomialVector, lib::PolyhedraLibrary)
    p = newtonpolytope(Z, lib)
    mindeg, maxdeg = cfld(extdeg(Z), 2)
    MonomialVector(vars(Z), mindeg:maxdeg, z -> z in p)
end
getmonomialsforcertificate(Z::Vector, libs...) = getmonomialsforcertificate(MonomialVector(Z), libs...)
getmonomialsforcertificate(p::Polynomial, libs...) = getmonomialsforcertificate(monomials(p), libs...)

getmonomialsforcertificate(porZ) = getmonomialsforcertificate(porZ, CheapOuterLibrary())

function randpsd(n; r=n, eps=0.1)
    Q = randn(n,n)
    d = zeros(Float64, n)
    d[1:r] = eps + abs.(randn(r))
    Q' * Diagonal(d) * Q
end

function randsos(Z::MonomialVector; r=-1, monotype=:Classic, eps=0.1)
    if monotype == :Classic
        x = getmonomialsforcertificate(Z)
    elseif monotype == :Gram
        x = Z
    else
        throw(ArgumentError("Monotype $monotype not known"))
    end
    n = length(x)
    if r < 0
        r = n
    end
    MatPolynomial(randpsd(n, r=r, eps=eps), x)
end

randsos(Z::Vector; r=-1, monotype=:Classic, eps=0.1) = randsos(MonomialVector(Z), r=r, monotype=monotype, eps=eps)
