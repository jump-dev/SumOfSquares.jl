export SOSDecomposition, SOSDecompositionWithDomain, sos_decomposition

"""
    struct SOSDecomposition{T, PT}

Represents a Sum-of-Squares decomposition without domain.
"""
struct SOSDecomposition{T,PT<:_APL{T},U} <: AbstractDecomposition{U}
    ps::Vector{PT}
    function SOSDecomposition{T,PT,U}(ps::Vector{PT}) where {T,PT,U}
        return new(ps)
    end
end

function SOSDecomposition(ps::Vector{PT}) where {T,PT<:_APL{T}}
    return SOSDecomposition{T,PT,_promote_add_mul(T)}(ps)
end
function MP.polynomial_type(
    ::Union{SOSDecomposition{T,PT,U},Type{SOSDecomposition{T,PT,U}}},
) where {T,PT,U}
    return MP.polynomial_type(PT, U)
end

#function SOSDecomposition(ps::Vector)
#    T = reduce(promote_type, Int, map(eltype, ps))
#    SOSDecomposition{T}(ps)
#end

function GramMatrix(p::SOSDecomposition{T}) where {T}
    X = MP.merge_monomial_vectors(map(MP.monomials, p))
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
    return GramMatrix(Q' * Q, X)
end

function SOSDecomposition(
    p::GramMatrix,
    ranktol = 0.0,
    dec::MultivariateMoments.LowRankLDLTAlgorithm = SVDLDLT(),
)
    n = length(p.basis)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    ldlt =
        MultivariateMoments.low_rank_ldlt(Matrix(value_matrix(p)), dec, ranktol)
    ps = [
        MP.polynomial(√ldlt.singular_values[i] * ldlt.L[:, i], p.basis) for
        i in axes(ldlt.L, 2)
    ]
    return SOSDecomposition(ps)
end
# Without LDL^T, we need to do float(T)
#SOSDecomposition(p::GramMatrix{C, T}) where {C, T} = SOSDecomposition{C, float(T)}(p)

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
        MultivariateMoments.compare_modulo_permutation(
            (i, j) -> isapprox(p.ps[i], q.ps[j]; kwargs...),
            m,
        )
    end
end

function Base.promote_rule(
    ::Type{SOSDecomposition{T1,PT1,U1}},
    ::Type{SOSDecomposition{T2,PT2,U2}},
) where {T1,T2,PT1<:_APL{T1},PT2<:_APL{T2},U1,U2}
    T = promote_type(T1, T2)
    return SOSDecomposition{T,promote_type(PT1, PT2),_promote_add_mul(T)}
end

function Base.convert(
    ::Type{SOSDecomposition{T,PT,U}},
    p::SOSDecomposition,
) where {T,PT,U}
    return SOSDecomposition(convert(Vector{PT}, p.ps))
end

function MP.polynomial(decomp::SOSDecomposition)
    return sum(decomp.ps .^ 2)
end
function MP.polynomial(decomp::SOSDecomposition, T::Type)
    return MP.polynomial(MP.polynomial(decomp), T)
end

"""
    struct SOSDecompositionWithDomain{T, PT, S}

Represents a Sum-of-Squares decomposition on a basic semi-algebraic domain.
"""
struct SOSDecompositionWithDomain{T,PT<:_APL{T},U,S<:AbstractSemialgebraicSet}
    sos::SOSDecomposition{T,PT,U}
    sosj::Vector{SOSDecomposition{T,PT,U}}
    domain::S
end

function SOSDecompositionWithDomain(
    ps::SOSDecomposition{T1,PT1,U1},
    vps::Vector{SOSDecomposition{T2,PT2,U2}},
    set::AbstractSemialgebraicSet,
) where {T1,T2,PT1,PT2,U1,U2}
    ptype =
        promote_type(SOSDecomposition{T1,PT1,U1}, SOSDecomposition{T2,PT2,U2})
    return SOSDecompositionWithDomain(
        convert(ptype, ps),
        convert(Vector{ptype}, vps),
        set,
    )
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
        p += MP.polynomial(Gj) * gj
    end
    return p
end

function Base.isapprox(
    p::SOSDecompositionWithDomain,
    q::SOSDecompositionWithDomain;
    kwargs...,
)
    return isapprox(p.sos, q.sos) &&
           all(isapprox.(p.sosj, q.sosj)) &&
           p.domain == q.domain
end
