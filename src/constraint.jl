export certificate_basis, certificate_monomials, gram_matrix, lagrangian_multipliers
export NonnegPolyInnerCone, DSOSCone, SDSOSCone, SOSCone
export PSDMatrixInnerCone, SOSMatrixCone
export ConvexPolyInnerCone, SOSConvexCone

"""
    struct NonnegPolyInnerCone{MCT <: MOI.AbstractVectorSet}
    end

Inner approximation of the cone of nonnegative polynomials defined by the set
of polynomials of the form
```julia
X^\\top Q X
```
where `X` is a vector of monomials and `Q` is a symmetric matrix that belongs to
the cone `psd_inner`.
"""
struct NonnegPolyInnerCone{MCT <: MOI.AbstractVectorSet} <: SOSLikeCone
end
matrix_cone_type(::Type{NonnegPolyInnerCone{MCT}}) where {MCT} = MCT

_wrap(::MIME"text/plain", s) = s
_wrap(::MIME"text/latex", s) = "\\text{ " * s * "}"

"""
    const SOSCone = NonnegPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}

Sum-of-squares cone; see [`NonnegPolyInnerCone`](@ref).
"""
const SOSCone = NonnegPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}
function JuMP.in_set_string(print_mode, ::SOSCone)
    return _wrap(print_mode, "is SOS")
end

"""
    const SDSOSCone = NonnegPolyInnerCone{ScaledDiagonallyDominantConeTriangle}

Scaled-diagonally-dominant-sum-of-squares cone; see [Definition 2, AM17] and
[`NonnegPolyInnerCone`](@ref).

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.
"""
const SDSOSCone = NonnegPolyInnerCone{ScaledDiagonallyDominantConeTriangle}
function JuMP.in_set_string(print_mode, ::SDSOSCone)
    return _wrap(print_mode, "is SDSOS")
end

"""
    const DSOSCone = NonnegPolyInnerCone{DiagonallyDominantConeTriangle}

Diagonally-dominant-sum-of-squares cone; see [Definition 2, AM17] and
[`NonnegPolyInnerCone`](@ref).

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.
"""
const DSOSCone = NonnegPolyInnerCone{DiagonallyDominantConeTriangle}
function JuMP.in_set_string(print_mode, ::DSOSCone)
    return _wrap(print_mode, "is DSOS")
end

function JuMP.reshape_set(set::SOSPolynomialSet, ::PolyJuMP.PolynomialShape)
    return Certificate.cone(set.certificate)
end

function default_ideal_certificate(domain, basis, cone, maxdegree, newton_polytope)
    if maxdegree === nothing
        return Certificate.Newton(cone, basis, newton_polytope)
    else
        return Certificate.MaxDegree(cone, basis, maxdegree)
    end
end
function default_ideal_certificate(domain::FixedVariablesSet, basis, cone, maxdegree, newton_polytope)
    return Certificate.Remainder(Certificate.Newton(cone, basis, newton_polytope))
end
function default_ideal_certificate(domain::FullSpace, basis, cone, maxdegree, newton_polytope)
    return Certificate.Newton(cone, basis, newton_polytope)
end

function default_ideal_certificate(
    domain::AbstractAlgebraicSet, sparsity::Certificate.Sparsity.NoPattern, basis::AbstractPolynomialBasis, cone, args...)
    return Certificate.FixedBasis(cone, basis)
end
function default_ideal_certificate(
    domain::AbstractAlgebraicSet, sparsity::Certificate.Sparsity.NoPattern, args...)
    return default_ideal_certificate(domain, args...)
end
function default_ideal_certificate(
    domain::AbstractAlgebraicSet, sparsity::Certificate.Sparsity.Pattern, args...)
    return Certificate.Sparsity.Ideal(sparsity, args...)
end

function default_ideal_certificate(
    domain::AbstractAlgebraicSet, symmetry::Nothing, args...)
    return default_ideal_certificate(domain, args...)
end
function default_ideal_certificate(
    domain::AbstractAlgebraicSet, symmetry::Certificate.Symmetry.Pattern, args...)
    return Certificate.Symmetry.Ideal(symmetry, default_ideal_certificate(domain, args...))
end

function default_ideal_certificate(
    domain::AbstractAlgebraicSet, newton_of_remainder::Bool, args...)
    c = default_ideal_certificate(domain, args...)
    if newton_of_remainder && !(c isa SumOfSquares.Certificate.Remainder)
        return SumOfSquares.Certificate.Remainder(c)
    end
    return c
end

function default_ideal_certificate(domain::BasicSemialgebraicSet, args...)
    return default_ideal_certificate(domain.V, args...)
end

function default_certificate(::AbstractAlgebraicSet, sparsity, ideal_certificate, cone, basis, maxdegree)
    return ideal_certificate
end
function default_certificate(::BasicSemialgebraicSet, sparsity::Certificate.Sparsity.Pattern,
                             ideal_certificate::Certificate.Sparsity.Ideal, cone,
                             basis, maxdegree)
    return Certificate.Sparsity.Preorder(
        sparsity, Certificate.Putinar(ideal_certificate.certificate, cone, basis, maxdegree))
end
function default_certificate(::BasicSemialgebraicSet, ::Certificate.Sparsity.NoPattern,
                             ideal_certificate, cone, basis, maxdegree)
    return Certificate.Putinar(
        ideal_certificate, cone, basis, maxdegree)
end

# Julia v1.0 does not support `init` keyword
#_max_maxdegree(p) = maximum(MP.maxdegree, p, init=0)
_max_maxdegree(p) = mapreduce(MP.maxdegree, max, p, init=0)

_maxdegree(domain) = 0

function _maxdegree(domain::AlgebraicSet)
    return _max_maxdegree(domain.I.p)
end

function _maxdegree(domain::BasicSemialgebraicSet)
    return max(_max_maxdegree(domain.p), _maxdegree(domain.V))
end

"""
    default_maxdegree(monos, domain)

Return the default `maxdegree` to use for certifying a polynomial with
monomials `monos` to be Sum-of-Squares over the domain `domain`.
By default, the maximum of the maxdegree of `monos` and of all multipliers
of the domain are used so that at least constant multipliers can be used
with a Putinar certificate.
"""
function default_maxdegree(monos, domain)
    return max(MP.maxdegree(monos), _maxdegree(domain))
end

function JuMP.moi_set(
    cone::SOSLikeCone,
    monos::AbstractVector{<:MP.AbstractMonomial};
    domain::AbstractSemialgebraicSet=FullSpace(),
    basis=MonomialBasis,
    newton_polytope::Tuple=tuple(),
    maxdegree::Union{Nothing, Int}=default_maxdegree(monos, domain),
    sparsity::Certificate.Sparsity.Pattern=Certificate.Sparsity.NoPattern(),
    symmetry::Union{Nothing,Certificate.Symmetry.Pattern}=nothing,
    newton_of_remainder::Bool=false,
    ideal_certificate=default_ideal_certificate(
        domain, newton_of_remainder, symmetry, sparsity, basis, cone, maxdegree, newton_polytope),
    certificate=default_certificate(
        domain, sparsity, ideal_certificate, cone, basis, maxdegree)
    )
    return SOSPolynomialSet(domain, monos, certificate)
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{EmptyCone})
    return [Bridges.Constraint.EmptyBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{PositiveSemidefinite2x2ConeTriangle})
    return [Bridges.Constraint.PositiveSemidefinite2x2Bridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:DiagonallyDominantConeTriangle})
    return [Bridges.Constraint.DiagonallyDominantBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:ScaledDiagonallyDominantConeTriangle})
    return [Bridges.Constraint.ScaledDiagonallyDominantBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:SOSPolynomialSet{<:AbstractAlgebraicSet}})
    return [Bridges.Constraint.SOSPolynomialBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:SOSPolynomialSet{<:BasicSemialgebraicSet}})
    return [Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge]
end

function JuMP.build_constraint(_error::Function, p, cone::SOSLikeCone; kws...)
    coefs = PolyJuMP.non_constant_coefficients(p)
    monos = MP.monomials(p)
    set = JuMP.moi_set(cone, monos; kws...)
    shape = PolyJuMP.PolynomialShape(monos)
    return PolyJuMP.bridgeable(JuMP.VectorConstraint(coefs, set, shape),
                               JuMP.moi_function_type(typeof(coefs)),
                               typeof(set))
end

_non_constant(a::Vector{T}) where T = convert.(MOI.ScalarAffineFunction{T}, a)
_non_constant(a::Vector{<:MOI.AbstractFunction}) = a

# Add constraint with `p` having coefficients being MOI functions.
# This is needed as a workaround as JuMP does not support complex numbers yet.
# We can remove it when https://github.com/jump-dev/JuMP.jl/pull/2391 is done
# We overload `JuMP.add_constraint` to avoid clash with the name.
function JuMP.add_constraint(model::MOI.ModelLike, p, cone::SOSLikeCone; kws...)
    coefs = MOI.Utilities.vectorize(_non_constant(MP.coefficients(p)))
    monos = MP.monomials(p)
    set = JuMP.moi_set(cone, monos; kws...)
    return MOI.add_constraint(model, coefs, set)
end
function JuMP.add_constraint(model::JuMP.Model, p, cone::SOSLikeCone; kws...)
    ci = JuMP.add_constraint(JuMP.backend(model), p, cone; kws...)
    return JuMP.ConstraintRef(model, ci, PolyJuMP.PolynomialShape(MP.monomials(p)))
end

struct ValueNotSupported <: Exception end
function Base.showerror(io::IO, ::ValueNotSupported)
    print(io, "`value` is no supported for Sum-of-Squares constraints, use",
          " `gram_matrix` instead.")
end

struct DualNotSupported <: Exception end
function Base.showerror(io::IO, ::DualNotSupported)
    print(io, "`dual` is no supported for Sum-of-Squares constraints in a",
          " domain, use `moment_matrix` instead.")
end

"""
    gram_matrix(cref::JuMP.ConstraintRef)

Return the [`GramMatrixAttribute`](@ref) of `cref`.
"""
function gram_matrix(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, GramMatrixAttribute(), cref)
end

"""
    sos_decomposition(cref::JuMP.ConstraintRef)

Return the [`SOSDecompositionAttribute`](@ref) of `cref`.
"""
function sos_decomposition(cref::JuMP.ConstraintRef, ranktol::Real=0.0,
                           dec::MultivariateMoments.LowRankChol=SVDChol())
    return MOI.get(cref.model, SOSDecompositionAttribute(ranktol, dec), cref)
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

"""
    moment_matrix(cref::JuMP.ConstraintRef)

Return the [`MomentMatrixAttribute`](@ref) of `cref`.
"""
function MultivariateMoments.moment_matrix(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, MomentMatrixAttribute(), cref)
end

# Equivalent but more efficient than `moment_matrix(cref).basis` as it does not
# need to query any result from the solver
"""
    certificate_basis(cref::JuMP.ConstraintRef)

Return the [`CertificateBasis`](@ref) of `cref`.
"""
function certificate_basis(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, CertificateBasis(), cref)
end

"""
    certificate_monomials(cref::JuMP.ConstraintRef)

Return the monomials of [`certificate_basis`](@ref). If the basis if not
`MultivariateBases.AbstractMonomialBasis`, an error is thrown.
"""
function certificate_monomials(cref::JuMP.ConstraintRef)
    return basis_monomials(certificate_basis(cref))
end
basis_monomials(basis::AbstractMonomialBasis) = basis.monomials
basis_monomials(basis::AbstractPolynomialBasis) = error("`certificate_monomials` is not supported with `$(typeof(basis))`, use `certificate_basis` instead.")

"""
    lagrangian_multipliers(cref::JuMP.ConstraintRef)

Return the [`LagrangianMultipliers`](@ref) of `cref`.
"""
function lagrangian_multipliers(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, LagrangianMultipliers(), cref)
end

"""
    struct PSDMatrixInnerCone{MCT <: MOI.AbstractVectorSet} <: PolyJuMP.PolynomialSet
    end

Inner approximation of the cone of polynomial matrices `P(x)` that are positive
semidefinite for any `x` defined by the set of ``n \\times n`` polynomial
matrices such that the polynomial ``y^\\top P(x) y`` belongs to
`NonnegPolyInnerCone{MCT}` where `y` is a vector of ``n`` auxiliary polynomial
variables.
"""
struct PSDMatrixInnerCone{MCT <: MOI.AbstractVectorSet} <: PolyJuMP.PolynomialSet
end

"""
    const SOSMatrixCone = PSDMatrixInnerCone{MOI.PositiveSemidefiniteConeTriangle}

Sum-of-squares matrices cone; see [Section 3.3.2, BPT12] and
[`PSDMatrixInnerCone`](@ref).

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, 2012.
"""
const SOSMatrixCone = PSDMatrixInnerCone{MOI.PositiveSemidefiniteConeTriangle}

function JuMP.build_constraint(_error::Function, P::AbstractMatrix{<:MP.APL},
                               ::PSDMatrixInnerCone{MCT};
                               newton_polytope::Tuple = tuple(),
                               kws...) where MCT
    n = LinearAlgebra.checksquare(P)
    if !issymmetric(P)
        _error("The polynomial matrix constrained to be SOS must be symmetric.")
    end
    y = [MP.similarvariable(eltype(P), gensym()) for i in 1:n]
    p = dot(y, P * y)
    # TODO Newton_polytope=(y,) may not be the best idea if exact newton
    #      polytope computation is used.
    #      See "Sum-of-Squares Matrices" notebook
    JuMP.build_constraint(_error, p, NonnegPolyInnerCone{MCT}();
                          newton_polytope=(y, newton_polytope...), kws...)
end

"""
    struct ConvexPolyInnerCone{MCT} <: PolyJuMP.PolynomialSet end

Inner approximation of the set of convex polynomials defined by the set of
polynomials such that their hessian belongs to to the set
`PSDMatrixInnerCone{MCT}()`
"""
struct ConvexPolyInnerCone{MCT} <: PolyJuMP.PolynomialSet end

"""
    const SOSConvexCone = ConvexPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}

Sum-of-squares convex polynomials cone; see [Section 3.3.3, BPT12] and
[`ConvexPolyInnerCone`](@ref).

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, 2012.
"""
const SOSConvexCone = ConvexPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}

function JuMP.build_constraint(_error::Function, p::MP.APL,
                               ::ConvexPolyInnerCone{MCT}; kws...) where MCT
    hessian = MP.differentiate(p, MP.variables(p), 2)
    JuMP.build_constraint(_error, hessian, PSDMatrixInnerCone{MCT}(); kws...)
end
