using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl#bl/nonperm"

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using DynamicPolynomials
@polyvar x

poly = x^4 - 2x^2

using SymbolicWedderburn
using PermutationGroups
using Cyclotomics
using SumOfSquares

struct ScaledPerm{T, I} <: AbstractPerm
    indices::Vector{Pair{I, T}}
end
function trace_matrix_representative(p::ScaledPerm)
    return sum(pair.second for (i, pair) in enumerate(p.indices) if pair.first == i)
end
function Base.inv(p::ScaledPerm{T,I}) where {T,I}
    indices = Vector{Pair{I,T}}(undef, length(p.indices))
    for i in eachindex(indices)
        pair = p.indices[i]
        indices[pair.first] = i => inv(pair.second)
    end
    return ScaledPerm(indices)
end

function PermutationGroups.order(p::ScaledPerm)
    cur = p
    _one = one(p)
    o = 1
    while cur != _one
        o += 1
        cur *= p
    end
    return o
end
Base.one(p::ScaledPerm{T}) where {T} = ScaledPerm([i => one(T) for i in eachindex(p.indices)])
Base.:(==)(p::ScaledPerm, q::ScaledPerm) = p.indices == q.indices
Base.hash(p::ScaledPerm, u::UInt64) = hash(p.indices, u)
SymbolicWedderburn.degree(p::ScaledPerm) = length(p.indices)
Base.:^(i::Integer, p::ScaledPerm) = p.indices[i]
function SymbolicWedderburn.add_inverse_permutation!(result, val, i::Int, j::Pair)
    result[i, j.first] += val / j.second
end
function SymbolicWedderburn.action_character(conjugacy_cls::AbstractVector{<:AbstractOrbit{<:ScaledPerm}})
    vals = Int[trace_matrix_representative(first(cc)) for cc in conjugacy_cls]
    return SymbolicWedderburn.Character(vals, conjugacy_cls)
end
function Base.:*(p::ScaledPerm, q::ScaledPerm)
    return ScaledPerm(map(eachindex(q.indices)) do i
        pair_q = q.indices[i]
        pair_p = p.indices[pair_q.first]
        pair_p.first => pair_p.second * pair_q.second
    end)
end

function permutation(ehom::SymbolicWedderburn.ExtensionHomomorphism{T}, els::Vector{T}) where {T}
    return Perm([ehom[el] for el in els])
end
function permutation(ehom::SymbolicWedderburn.ExtensionHomomorphism{<:AbstractMonomial}, terms::Vector{<:AbstractTerm})
    return ScaledPerm([ehom[monomial(term)] => coefficient(term) for term in terms])
end
function (ehom::SymbolicWedderburn.ExtensionHomomorphism)(action)
    return permutation(ehom, [ehom.op(f, action) for f in ehom.features])
end

function SymbolicWedderburn.ExtensionHomomorphism(basis::MB.MonomialBasis, action)
    monos = collect(basis.monomials)
    mono_to_index = Dict(monos[i] => i for i in eachindex(monos))
    return SymbolicWedderburn.ExtensionHomomorphism(monos, mono_to_index, action)
end

function MP.polynomialtype(::Type{<:MB.AbstractPolynomialVectorBasis{PT}}, T::Type) where PT
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end

struct SymmetricIdeal{CT, GT} <: Certificate.AbstractIdealCertificate
    cone::CT
    group::GT
end
SumOfSquares.matrix_cone_type(::Type{<:SymmetricIdeal{CT}}) where {CT} = SumOfSquares.matrix_cone_type(CT)
Certificate.get(::Type{<:SymmetricIdeal}, ::SumOfSquares.Certificate.GramBasisType) = Vector{MB.FixedPolynomialBasis}
Certificate.zero_basis_type(::Type{<:SymmetricIdeal}) = MB.MonomialBasis
Certificate.zero_basis(::SymmetricIdeal) = MB.MonomialBasis
Certificate.get(::SymmetricIdeal, ::Certificate.ReducedPolynomial, poly, domain) = poly
function Certificate.get(cert::SymmetricIdeal, ::Certificate.GramBasis, poly)
    basis = Certificate.maxdegree_gram_basis(MB.MonomialBasis, MP.variables(poly), MP.maxdegree(poly))
    R = SymbolicWedderburn.symmetry_adapted_basis(Float64, cert.group, basis, action)
    return map(R) do Ri
        FixedPolynomialBasis(convert(Matrix{Float64}, Ri) * basis.monomials)
    end
end

struct EvenOddAction <: GroupElem
    identity::Bool
end
PermutationGroups.order(a::EvenOddAction) = a.identity ? 1 : 2
Base.one(::EvenOddAction) = EvenOddAction(true)
Base.inv(a::EvenOddAction) = a
PermutationGroups.mul!(::EvenOddAction, a::EvenOddAction, b::EvenOddAction) = EvenOddAction(xor(a.identity, b.identity))
function action(mono::AbstractMonomial, a::EvenOddAction)
    if a.identity || iseven(MP.degree(mono))
        return 1 * mono
    else
        return -1 * mono
    end
end

struct EvenOddSymmetry <: Group
end
_orbit(cc) = PermutationGroups.Orbit(cc, Dict(a => nothing for a in cc))
SymbolicWedderburn.conjugacy_classes_orbit(::EvenOddSymmetry) = [_orbit([EvenOddAction(true)]), _orbit([EvenOddAction(false)])]

G = EvenOddSymmetry()

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G))
optimize!(model)
@show value(t)
value(t)

gram_matrix(con_ref)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

