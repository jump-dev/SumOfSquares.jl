export SOSDecomposition, SOSDecompositionWithDomain, sos_decomposition

"""
    struct SOSDecomposition{T, PT}

Represents a Sum-of-Squares decomposition without domain.
"""
struct SOSDecomposition{A,T,V,U} <: AbstractDecomposition{U}
    ps::Vector{SA.AlgebraElement{A,T,V}} # TODO rename `elements`
    function SOSDecomposition{A,T,V,U}(
        ps::Vector{SA.AlgebraElement{A,T,V}},
    ) where {A,T,V,U}
        return new(ps)
    end
end

function SOSDecomposition(
    elements::Vector{SA.AlgebraElement{A,T,V}},
) where {A,T,V}
    return SOSDecomposition{A,T,V,_promote_add_mul(T)}(elements)
end
function MP.polynomial_type(
    ::Union{SOSDecomposition{A,T,V,U},Type{SOSDecomposition{A,T,V,U}}},
) where {A,T,V,U}
    return MP.polynomial_type(MP.polynomial_type(SA.AlgebraElement{A,T,V}), U)
end

function GramMatrix(p::SOSDecomposition{A,T}) where {A,T}
    basis = mapreduce(SA.basis, (b1, b2) -> MB.merge_bases(b1, b2)[1], p.ps)
    m = length(p.ps)
    n = length(basis)
    Q = zeros(T, m, n)
    for i in 1:m
        j = 1
        for (k, v) in SA.nonzero_pairs(SA.coeffs(p.ps[i]))
            poly = SA.basis(p.ps[i])[k]
            while j in eachindex(basis) && basis[j] != poly
                j += 1
            end
            Q[i, j] = v
            j += 1
        end
    end
    return GramMatrix(Q' * Q, basis)
end

_lazy_adjoint(x::AbstractVector{<:Real}) = x
_lazy_adjoint(x::AbstractVector) = adjoint.(x)

function SOSDecomposition(
    p::GramMatrix,
    ranktol = 0.0,
    dec::MultivariateMoments.LowRankLDLTAlgorithm = SVDLDLT(),
)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    ldlt =
        MultivariateMoments.low_rank_ldlt(Matrix(value_matrix(p)), dec, ranktol)
    # The Sum-of-Squares decomposition is
    # ∑ adjoint(u_i) * u_i
    # and we have `L` of the LDL* so we need to take the adjoint.
    return SOSDecomposition(
        map(axes(ldlt.L, 2)) do i
            return MB.algebra_element(
                √ldlt.singular_values[i] * _lazy_adjoint(ldlt.L[:, i]),
                p.basis,
            )
        end,
    )
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
