export certificate_basis,
    certificate_monomials, gram_matrix, lagrangian_multipliers
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
struct NonnegPolyInnerCone{MCT<:MOI.AbstractVectorSet} <: SOSLikeCone end
matrix_cone_type(::Type{NonnegPolyInnerCone{MCT}}) where {MCT} = MCT

_wrap(::MIME"text/plain", s) = s
_wrap(::MIME"text/latex", s) = "\\text{ " * s * "}"

"""
    const SOSCone = NonnegPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}

Sum-of-squares cone; see [`NonnegPolyInnerCone`](@ref).
"""
const SOSCone = NonnegPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}
function JuMP.in_set_string(print_mode::MIME, ::SOSCone)
    return _wrap(print_mode, "is SOS")
end

"""
    const SDSOSCone = NonnegPolyInnerCone{ScaledDiagonallyDominantConeTriangle}

Scaled-diagonally-dominant-sum-of-squares cone; see [Ahmadi2017; Definition 2](@cite) and
[`NonnegPolyInnerCone`](@ref).
"""
const SDSOSCone = NonnegPolyInnerCone{ScaledDiagonallyDominantConeTriangle}
function JuMP.in_set_string(print_mode::MIME, ::SDSOSCone)
    return _wrap(print_mode, "is SDSOS")
end

"""
    const DSOSCone = NonnegPolyInnerCone{DiagonallyDominantConeTriangle}

Diagonally-dominant-sum-of-squares cone; see [Ahmadi2017; Definition 2](@cite) and
[`NonnegPolyInnerCone`](@ref).
"""
const DSOSCone = NonnegPolyInnerCone{DiagonallyDominantConeTriangle}
function JuMP.in_set_string(print_mode::MIME, ::DSOSCone)
    return _wrap(print_mode, "is DSOS")
end

function JuMP.reshape_set(set::SOSPolynomialSet, ::PolyJuMP.PolynomialShape)
    return Certificate.cone(set.certificate)
end

function default_ideal_certificate(
    ::AbstractAlgebraicSet,
    gram_basis::SA.ImplicitBasis,
    zero_basis,
    cone,
    maxdegree,
    newton_polytope,
)
    if maxdegree === nothing
        return Certificate.Newton(cone, gram_basis, zero_basis, newton_polytope)
    else
        return Certificate.MaxDegree(cone, gram_basis, zero_basis, maxdegree)
    end
end
function default_ideal_certificate(
    ::FixedVariablesSet,
    gram_basis::SA.ImplicitBasis,
    zero_basis,
    cone,
    maxdegree,
    newton_polytope,
)
    return Certificate.Remainder(
        Certificate.Newton(cone, gram_basis, zero_basis, newton_polytope),
    )
end
function default_ideal_certificate(
    ::FullSpace,
    gram_basis::SA.ImplicitBasis,
    zero_basis,
    cone,
    maxdegree,
    newton_polytope,
)
    if newton_polytope === nothing
        return Certificate.MaxDegree(cone, gram_basis, zero_basis, maxdegree)
    else
        return Certificate.Newton(cone, gram_basis, zero_basis, newton_polytope)
    end
end

function default_ideal_certificate(
    ::AbstractAlgebraicSet,
    gram_basis::SA.ExplicitBasis,
    zero_basis,
    cone,
    args...,
)
    return Certificate.FixedBasis(cone, gram_basis, zero_basis)
end
function default_ideal_certificate(
    domain::AbstractAlgebraicSet,
    ::Certificate.Sparsity.NoPattern,
    args...,
)
    return default_ideal_certificate(domain, args...)
end
function default_ideal_certificate(
    ::AbstractAlgebraicSet,
    sparsity::Certificate.Sparsity.Pattern,
    args...,
)
    return Certificate.Sparsity.Ideal(sparsity, args...)
end

function default_ideal_certificate(
    domain::AbstractAlgebraicSet,
    ::Nothing,
    args...,
)
    return default_ideal_certificate(domain, args...)
end
function default_ideal_certificate(
    domain::AbstractAlgebraicSet,
    symmetry::Certificate.Symmetry.Pattern,
    args...,
)
    return Certificate.Symmetry.Ideal(
        symmetry,
        default_ideal_certificate(domain, args...),
    )
end

function default_ideal_certificate(
    domain::AbstractAlgebraicSet,
    newton_of_remainder::Bool,
    args...,
)
    c = default_ideal_certificate(domain, args...)
    if newton_of_remainder && !(c isa SumOfSquares.Certificate.Remainder)
        return SumOfSquares.Certificate.Remainder(c)
    end
    return c
end

function default_ideal_certificate(domain::BasicSemialgebraicSet, args...)
    return default_ideal_certificate(
        SemialgebraicSets.algebraic_set(domain),
        args...,
    )
end

function default_certificate(
    ::AbstractAlgebraicSet,
    sparsity,
    ideal_certificate,
    cone,
    gram_basis,
    zero_basis,
    maxdegree,
    newton_polytope,
)
    return ideal_certificate
end
function default_certificate(
    domain::BasicSemialgebraicSet,
    sparsity::Certificate.Sparsity.Pattern,
    ideal_certificate::Certificate.Sparsity.Ideal,
    cone,
    gram_basis,
    zero_basis,
    maxdegree,
    newton_polytope,
)
    nonsparse = default_certificate(
        domain,
        Certificate.Sparsity.NoPattern(),
        ideal_certificate.certificate,
        cone,
        gram_basis,
        zero_basis,
        maxdegree,
        newton_polytope,
    )
    return Certificate.Sparsity.Preorder(sparsity, nonsparse)
end
function default_certificate(
    domain::BasicSemialgebraicSet,
    ::Certificate.Sparsity.NoPattern,
    ideal_certificate,
    cone,
    gram_basis,
    zero_basis,
    maxdegree,
    newton_polytope,
)
    # We could take `multipliers_certificate = ideal_certificate` here but
    # that wouldn't work if `ideal_certificate` is `Remainder`,
    # `Sparseity.Ideal` or `Symmetry.Ideal`
    multipliers_certificate = default_ideal_certificate(
        SemialgebraicSets.algebraic_set(domain),
        gram_basis,
        zero_basis,
        cone,
        maxdegree,
        newton_polytope,
    )
    if multipliers_certificate isa Certificate.Remainder
        # TODO not supported yet so we drop the `Remainder` part
        multipliers_certificate = multipliers_certificate.gram_certificate
    end
    return Certificate.Putinar(
        multipliers_certificate,
        ideal_certificate,
        maxdegree,
    )
end

# Julia v1.0 does not support `init` keyword
#_max_maxdegree(p) = maximum(MP.maxdegree, p, init=0)
_max_maxdegree(p) = mapreduce(MP.maxdegree, max, p, init = 0)

_maxdegree(domain) = 0

function _maxdegree(domain::AlgebraicSet)
    return _max_maxdegree(SemialgebraicSets.equalities(domain))
end

function _maxdegree(domain::BasicSemialgebraicSet)
    return max(
        _max_maxdegree(SemialgebraicSets.inequalities(domain)),
        _maxdegree(SemialgebraicSets.algebraic_set(domain)),
    )
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
    basis::SA.ExplicitBasis,
    gram_basis::SA.AbstractBasis,
    zero_basis::SA.AbstractBasis;
    domain::AbstractSemialgebraicSet = FullSpace(),
    newton_polytope::Union{Nothing,Tuple} = tuple(),
    maxdegree::Union{Nothing,Int} = default_maxdegree(basis, domain),
    sparsity::Certificate.Sparsity.Pattern = Certificate.Sparsity.NoPattern(),
    symmetry::Union{Nothing,Certificate.Symmetry.Pattern} = nothing,
    newton_of_remainder::Bool = false,
    ideal_certificate = default_ideal_certificate(
        domain,
        newton_of_remainder,
        symmetry,
        sparsity,
        gram_basis,
        zero_basis,
        cone,
        maxdegree,
        newton_polytope,
    ),
    certificate = default_certificate(
        domain,
        sparsity,
        ideal_certificate,
        cone,
        gram_basis,
        zero_basis,
        maxdegree,
        newton_polytope,
    ),
)
    return SOSPolynomialSet(domain, basis, certificate)
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{EmptyCone},
)
    return [(Bridges.Constraint.EmptyBridge, Float64)]
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{PositiveSemidefinite2x2ConeTriangle},
)
    return [(Bridges.Constraint.PositiveSemidefinite2x2Bridge, Float64)]
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:DiagonallyDominantConeTriangle},
)
    return [(Bridges.Constraint.DiagonallyDominantBridge, Float64)]
end

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ScaledDiagonallyDominantConeTriangle},
) # Needed so that `Variable.ScaledDiagonallyDominantBridge` is added as well
    return Tuple{Type,Type}[(
        MOI.Bridges.Constraint.VectorSlackBridge,
        PolyJuMP.coefficient_type_or_float(F),
    )]
end

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:WeightedSOSCone},
) # Needed so that `Variable.KernelBridge` is added as well
    return Tuple{Type,Type}[(
        MOI.Bridges.Constraint.VectorSlackBridge,
        PolyJuMP.coefficient_type_or_float(F),
    )]
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    S::Type{<:SOSPolynomialSet{<:AbstractAlgebraicSet}},
)
    return Tuple{Type,Type}[(
        Bridges.Constraint.SOSPolynomialBridge,
        _bridge_coefficient_type(S),
    )]
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    S::Type{<:SOSPolynomialSet{<:BasicSemialgebraicSet}},
)
    return [(
        Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge,
        _bridge_coefficient_type(S),
    )]
end

# Syntax: `@constraint(model, a >= b, SOSCone())`
function JuMP.build_constraint(
    _error::Function,
    f,
    ::JuMP.Nonnegatives,
    extra::SOSLikeCone;
    kws...,
)
    return build_constraint(_error, f, extra; kws...)
end

_promote_coef_type(::Type{V}, ::Type) where {V<:JuMP.AbstractVariableRef} = V
_promote_coef_type(::Type{F}, ::Type{T}) where {F,T} = promote_type(F, T)

function _default_basis(
    coeffs,
    basis::MB.SubBasis{B},
    gram_basis::MB.MonomialIndexedBasis{G},
) where {B,G}
    if B === G
        return coeffs, basis
    else
        return _default_basis(
            SA.SparseCoefficients(basis.monomials, coeffs),
            MB.implicit_basis(basis),
            gram_basis,
        )
    end
end

function _default_basis(
    p::SA.AbstractCoefficients,
    basis::MB.FullBasis{B},
    gram_basis::MB.MonomialIndexedBasis{G},
) where {B,G}
    if B === G
        return _default_basis(
            collect(SA.values(p)),
            MB.SubBasis{B}(collect(SA.keys(p))),
            gram_basis,
        )
    else
        new_basis = MB.FullBasis{G,MP.monomial_type(typeof(basis))}()
        return _default_basis(
            SA.coeffs(p, basis, new_basis),
            new_basis,
            gram_basis,
        )
    end
end

function _default_gram_basis(
    ::MB.MonomialIndexedBasis{B,M},
    ::Nothing,
) where {B,M}
    return MB.FullBasis{B,M}()
end

function _default_gram_basis(
    ::MB.MonomialIndexedBasis{_B,M},
    ::Type{B},
) where {_B,B,M}
    return MB.FullBasis{B,M}()
end

function _default_gram_basis(_, basis::MB.MonomialIndexedBasis)
    return basis
end

function _default_basis(a::SA.AlgebraElement, basis)
    b = SA.basis(a)
    gram = _default_gram_basis(b, basis)
    return _default_basis(SA.coeffs(a), b, gram)..., gram
end

function _default_basis(p::MP.AbstractPolynomialLike, basis)
    return _default_basis(
        MB.algebra_element(
            MB.sparse_coefficients(p),
            MB.FullBasis{MB.Monomial,MP.monomial_type(p)}(),
        ),
        basis,
    )
end

_default_zero_basis(basis, ::Nothing) = MB.implicit_basis(basis)
function _default_zero_basis(basis, nodes::MB.AbstractNodes)
    return MB.ImplicitLagrangeBasis(MP.variables(basis), nodes)
end

function JuMP.build_constraint(
    _error::Function,
    p,
    cone::SOSLikeCone;
    basis = nothing,
    zero_basis = nothing,
    kws...,
)
    __coefs, basis, gram_basis = _default_basis(p, basis)
    set = JuMP.moi_set(
        cone,
        basis,
        gram_basis,
        _default_zero_basis(basis, zero_basis);
        kws...,
    )
    _coefs = PolyJuMP.non_constant(__coefs)
    # If a polynomial with real coefficients is used with the Hermitian SOS
    # cone, we want to promote the coefficients to complex
    T = _bridge_coefficient_type(typeof(set))
    coefs = convert(Vector{_promote_coef_type(eltype(_coefs), T)}, _coefs)
    shape = PolyJuMP.PolynomialShape(basis)
    return PolyJuMP.bridgeable(
        JuMP.VectorConstraint(coefs, set, shape),
        JuMP.moi_function_type(typeof(coefs)),
        typeof(set),
    )
end

struct ValueNotSupported <: Exception end
function Base.showerror(io::IO, ::ValueNotSupported)
    return print(
        io,
        "`value` is no supported for Sum-of-Squares constraints, use",
        " `gram_matrix` instead.",
    )
end

struct DualNotSupported <: Exception end
function Base.showerror(io::IO, ::DualNotSupported)
    return print(
        io,
        "`dual` is no supported for Sum-of-Squares constraints in a",
        " domain, use `moment_matrix` instead.",
    )
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
function sos_decomposition(
    cref::JuMP.ConstraintRef,
    ranktol::Real = 0.0,
    dec::MultivariateMoments.LowRankLDLTAlgorithm = SVDLDLT(),
)
    return MOI.get(cref.model, SOSDecompositionAttribute(; ranktol, dec), cref)
end

"""
    sos_decomposition(cref::JuMP.ConstraintRef, K<:AbstractBasicSemialgebraicSet)

Return representation in the quadraic module associated with K.
"""
function sos_decomposition(
    cref::JuMP.ConstraintRef,
    K::AbstractBasicSemialgebraicSet,
    args...,
)
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
`MultivariateBases.SubBasis`, an error is thrown.
"""
function certificate_monomials(cref::JuMP.ConstraintRef)
    return basis_monomials(certificate_basis(cref))
end
basis_monomials(basis::MB.SubBasis) = basis.monomials
function basis_monomials(basis::SA.AbstractBasis)
    return error(
        "`certificate_monomials` is not supported with `$(typeof(basis))`, use `certificate_basis` instead.",
    )
end

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
struct PSDMatrixInnerCone{MCT<:MOI.AbstractVectorSet} <: PolyJuMP.PolynomialSet end

"""
    const SOSMatrixCone = PSDMatrixInnerCone{MOI.PositiveSemidefiniteConeTriangle}

Sum-of-squares matrices cone; see [Blekherman2012; Section 3.3.2](@cite) and
[`PSDMatrixInnerCone`](@ref).
"""
const SOSMatrixCone = PSDMatrixInnerCone{MOI.PositiveSemidefiniteConeTriangle}

function JuMP.build_constraint(
    _error::Function,
    P::AbstractMatrix{<:_AE_APL},
    ::PSDMatrixInnerCone{MCT};
    newton_polytope::Tuple = tuple(),
    kws...,
) where {MCT}
    n = LinearAlgebra.checksquare(P)
    if !issymmetric(P)
        _error("The polynomial matrix constrained to be SOS must be symmetric.")
    end
    y = [MP.similar_variable(eltype(P), gensym()) for i in 1:n]
    p = dot(y, P * y)
    # TODO Newton_polytope=(y,) may not be the best idea if exact newton
    #      polytope computation is used.
    #      See "Sum-of-Squares Matrices" notebook
    return JuMP.build_constraint(
        _error,
        p,
        NonnegPolyInnerCone{MCT}();
        newton_polytope = (y, newton_polytope...),
        kws...,
    )
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

Sum-of-squares convex polynomials cone; see [Blekherman2012; Section 3.3.3](@cite) and
[`ConvexPolyInnerCone`](@ref).
"""
const SOSConvexCone = ConvexPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}

function JuMP.build_constraint(
    _error::Function,
    p::_AE_APL,
    ::ConvexPolyInnerCone{MCT};
    kws...,
) where {MCT}
    hessian = MP.differentiate(p, MP.variables(p), 2)
    return JuMP.build_constraint(
        _error,
        hessian,
        PSDMatrixInnerCone{MCT}();
        kws...,
    )
end
