function SymbolicWedderburn.decompose(
    k::MP.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if features(hom) == monomials

    indcs = [hom[mono] for mono in MP.monomials(k)]
    coeffs = MP.coefficients(k)

    return indcs, coeffs
end

function SymbolicWedderburn.ExtensionHomomorphism(
    action::SymbolicWedderburn.Action,
    basis::MB.MonomialBasis,
)
    monos = collect(basis.monomials)
    mono_to_index = Dict(monos[i] => i for i in eachindex(monos))
    return SymbolicWedderburn.ExtensionHomomorphism(action, monos, mono_to_index)
end

struct VariablePermutation <: SymbolicWedderburn.ByPermutations end
_map_idx(f, v::AbstractVector) = map(f, eachindex(v))
_tuple_map_idx(f, ::Tuple{}, i) = tuple()
_tuple_map_idx(f, v::Tuple, i) = (f(i), _tuple_map_idx(f, Base.tail(v), i + 1)...)
_map_idx(f, v::Tuple) = _tuple_map_idx(f, v, 1)
function SymbolicWedderburn.action(::VariablePermutation, p, mono::MP.AbstractMonomial)
    v = MP.variables(mono)
    MP.substitute(MP.Eval(), mono, v => _map_idx(i -> v[i^p], v))
end
abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el,
    term::MP.AbstractTerm,
)
    return MP.coefficient(term) * SymbolicWedderburn.action(a, el, MP.monomial(term))
end
function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el,
    poly::MP.AbstractPolynomial,
)
    return sum([SymbolicWedderburn.action(a, el, term) for term in MP.terms(poly)])
end

# TODO Move it to MultivariateBases
function MP.polynomialtype(
    ::Type{<:MB.AbstractPolynomialVectorBasis{PT}},
    T::Type,
) where {PT}
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end

struct Ideal{C,GT,AT<:SymbolicWedderburn.Action} <: SumOfSquares.Certificate.AbstractIdealCertificate
    pattern::Pattern{GT,AT}
    certificate::C
end
SumOfSquares.Certificate.get(certificate::Ideal, attr::SumOfSquares.Certificate.Cone) = SumOfSquares.Certificate.get(certificate.certificate, attr)
function SumOfSquares.matrix_cone_type(::Type{<:Ideal{C}}) where {C}
    return SumOfSquares.matrix_cone_type(C)
end
SumOfSquares.Certificate.get(::Type{<:Ideal}, ::SumOfSquares.Certificate.GramBasisType) = Vector{Vector{MB.FixedPolynomialBasis}}
SumOfSquares.Certificate.zero_basis_type(::Type{<:Ideal}) = MB.MonomialBasis
SumOfSquares.Certificate.zero_basis(::Ideal) = MB.MonomialBasis
function SumOfSquares.Certificate.get(certificate::Ideal, attr::SumOfSquares.Certificate.ReducedPolynomial, poly, domain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, poly, domain)
end
function SumOfSquares.Certificate.get(cert::Ideal, attr::SumOfSquares.Certificate.GramBasis, poly)
    basis = SumOfSquares.Certificate.get(cert.certificate, attr, poly)
    T = SumOfSquares._complex(Float64, SumOfSquares.matrix_cone_type(typeof(cert)))
    summands = SymbolicWedderburn.symmetry_adapted_basis(T, cert.pattern.group, basis, cert.pattern.action)
    return map(summands) do summand
        R = SymbolicWedderburn.basis(summand)
        m = SymbolicWedderburn.multiplicity(summand)
        F = convert(Matrix{T}, R)
        N = size(R, 1)
        d = SymbolicWedderburn.degree(summand)
        if d > 1
            if m > 1
                ps = R * basis.monomials
                S = map(SymbolicWedderburn.gens(cert.pattern.group)) do g
                    Si = Matrix{T}(undef, N, N)
                    for i in eachindex(ps)
                        p = ps[i]
                        q = SymbolicWedderburn.action(cert.pattern.action, g, p)
                        coefs = MP.coefficients(q, basis.monomials)
                        col = row_echelon_linsolve(R, coefs)
                        Si[:, i] = col
                    end
                    return Si
                end
                U = ordered_block_diag(S, d)
            else
                U = Matrix{T}(LinearAlgebra.I, N, N)
            end
            map(1:d) do i
                MB.FixedPolynomialBasis(
                    (transpose(U[:, i:d:(i + d * (m - 1))]) * F) * basis.monomials,
                )
            end
        else
            [MB.FixedPolynomialBasis(F * basis.monomials)]
        end
    end
end
