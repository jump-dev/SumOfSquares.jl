export SOSDecomposition, SOSDecompositionWithDomain, sos_decomposition

"""
    struct SOSDecomposition{T, PT}

Represend SOSDecomposition without domain
"""
struct SOSDecomposition{T, PT <: MP.APL{T}} <: MP.APL{T} # If promote_op((x, y) -> x * y + x * y, T, T) != T then it might not be true
    ps::Vector{PT}
    function SOSDecomposition{T, PT}(ps::Vector{PT}) where {T, PT}
        new(ps)
    end
end

SOSDecomposition(ps::Vector{PT}) where {T, PT <: MP.APL{T}} = SOSDecomposition{T, PT}(ps)
MP.polynomialtype(::Type{SOSDecomposition{T, PT}}) where {T, PT} = MP.polynomialtype(PT)

#function SOSDecomposition(ps::Vector)
#    T = reduce(promote_type, Int, map(eltype, ps))
#    SOSDecomposition{T}(ps)
#end

function GramMatrix(p::SOSDecomposition{T}) where {T}
    X = MP.mergemonovec(map(MP.monomials, p))
    m = length(p)
    n = length(X)
    Q = zeros(T, m, n)
    for i in 1:m
        j = 1
        for t in MP.terms(p[i])
            while X[j] != MP.monomial(t)
                j += 1
            end
            Q[i, j] = MP.coefficient(t)
            j += 1
        end
    end
    GramMatrix(Q' * Q, X)
end

function SOSDecomposition(p::GramMatrix, ranktol=0.0,
                          dec::MultivariateMoments.LowRankChol=SVDChol())
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    nM, cM, Q = MultivariateMoments.lowrankchol(Matrix(getmat(p)), dec, ranktol)
    ps = [MP.polynomial(Q[i,:], p.x) for i in 1:size(Q, 1)]
    return SOSDecomposition(ps)
end
# Without LDL^T, we need to do float(T)
SOSDecomposition(p::GramMatrix{C, T}) where {C, T} = SOSDecomposition{C, float(T)}(p)

Base.length(p::SOSDecomposition) = length(p.ps)
Base.isempty(p::SOSDecomposition) = isempty(p.ps)
Base.iterate(p::SOSDecomposition, args...) = Base.iterate(p.ps, args...)
Base.getindex(p::SOSDecomposition, i::Int) = p.ps[i]

(p::GramMatrix)(s::MP.AbstractSubstitution...) = MP.polynomial(p)(s...)

function Base.show(io::IO, p::SOSDecomposition)
    for (i, q) in enumerate(p)
        print(io, "(")
        print(io, q)
        print(io, ")^2")
        if i != length(p)
            print(io, " + ")
        end
    end
end

function Base.isapprox(p::SOSDecomposition, q::SOSDecomposition; kwargs...)
    m = length(p.ps)
    if length(q.ps) != m
        false
    else
        MultivariateMoments.permcomp((i, j) -> isapprox(p.ps[i], q.ps[j]; kwargs...), m)
    end
end

function Base.promote_rule(::Type{SOSDecomposition{T1, PT1}}, ::Type{SOSDecomposition{T2, PT2}}) where {T1, T2, PT1<:MP.APL{T1}, PT2<:MP.APL{T2}}
    return SOSDecomposition{promote_type(T1, T2), promote_type(PT1, PT2)} 
end

function Base.convert(::Type{SOSDecomposition{T, PT}}, p::SOSDecomposition) where {T, PT}
    return SOSDecomposition(convert(Vector{PT}, p.ps))
end

function MP.polynomial(decomp::SOSDecomposition)
    return sum(decomp.ps.^2)
end
function MP.polynomial(decomp::SOSDecomposition, T::Type)
    return MP.polynomial(MP.polynomial(decomp), T)
end

"""
    function sos_decomposition(cref::JuMP.ConstraintRef, ranktol::Float64, dec::MultivariateMoments.lowrankchol)

Return representation as a sum of squares.
"""
function sos_decomposition(cref::JuMP.ConstraintRef, args...)
    return SOSDecomposition(gram_matrix(cref), args...)
end

"""
    struct SOSDecompositionWithDomain{T, PT, S}

Represend SOSDecomposition on a basic semi-algebraic domain.
"""
struct SOSDecompositionWithDomain{T, PT <: MP.APL{T}, S <: AbstractSemialgebraicSet }
    sos::SOSDecomposition{T, PT}
    sosj::Vector{SOSDecomposition{T, PT}}
    domain::S
end

function SOSDecompositionWithDomain(ps::SOSDecomposition{T1, PT1}, vps::Vector{SOSDecomposition{T2, PT2}}, set::AbstractSemialgebraicSet ) where {T1, T2, PT1, PT2}
    ptype = promote_type(SOSDecomposition{T1,PT1}, SOSDecomposition{T2, PT2})
    return SOSDecompositionWithDomain(convert(ptype, ps), convert(Vector{ptype}, vps), set)
end


function Base.show(io::IO, decomp::SOSDecompositionWithDomain)
    print(io, decomp.sos)
    for (sos, g) in zip(decomp.sosj, inequalities(decomp.domain))
        print(io, " + ")
        print(io, sos)
        print(io, " * ")
        print(io, "(")
        print(io, g)
        print(io, ")")
    end
end

function MP.polynomial(decomp::SOSDecompositionWithDomain)
    p = MP.polynomial(decomp.sos)
    if !(isempty(equalities(decomp.domain)))
        @error "Semialgebraic set has equality constraints"
    end
    for (Gj, gj) in zip(decomp.sosj, inequalities(decomp.domain))
        p += MP.polynomial(Gj)*gj
    end
    return p
end

function Base.isapprox(p::SOSDecompositionWithDomain, q::SOSDecompositionWithDomain; kwargs...)
    return isapprox(p.sos, q.sos) && all(isapprox.(p.sosj, q.sosj)) && p.domain == q.domain
end

"""
    sos_decomposition(cref::JuMP.ConstraintRef, K<:AbstractBasicSemialgebraicSet)

Return representation in the quadraic module associated with K. 
"""
function sos_decomposition(cref::JuMP.ConstraintRef, K::AbstractBasicSemialgebraicSet, args...)
    lm = SOSDecomposition.(lagrangian_multipliers(cref), args...)
    gm = sos_decomposition(cref, args...)
    return SOSDecompositionWithDomain(gm, lm, K)
end
