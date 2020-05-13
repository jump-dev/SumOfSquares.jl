export certificate_basis, certificate_monomials, gram_matrix, lagrangian_multipliers
export NonnegPolyInnerCone, DSOSCone, SDSOSCone, SOSCone
export PSDMatrixInnerCone, SOSMatrixCone
export ConvexPolyInnerCone, SOSConvexCone

export Sparsity, NoSparsity, VariableSparsity, MonomialSparsity
abstract type Sparsity end
struct NoSparsity <: Sparsity end
struct VariableSparsity <: Sparsity end
struct MonomialSparsity <: Sparsity end

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

_wrap(::Type{JuMP.REPLMode}, s) = s
_wrap(::Type{JuMP.IJuliaMode}, s) = "\\text{ " * s * "}"

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
    return Certificate.get(set.certificate, Certificate.Cone())
end
function default_ideal_certificate(domain, cone, basis, newton_polytope, maxdegree)
    if maxdegree === nothing
        return Certificate.Newton(cone, basis, newton_polytope)
    else
        return Certificate.MaxDegree(cone, basis, maxdegree)
    end
end
function default_ideal_certificate(domain::FixedVariablesSet, cone, basis, newton_polytope, maxdegree)
    return Certificate.Remainder(Certificate.Newton(cone, basis, newton_polytope))
end
function default_ideal_certificate(domain::FullSpace, cone, basis, newton_polytope, maxdegree)
    return Certificate.Newton(cone, basis, newton_polytope)
end

function default_ideal_certificate(
    domain::AbstractAlgebraicSet, cone, basis,
    newton_polytope, maxdegree, sparsity::Sparsity, newton_of_remainder::Bool)

    if sparsity isa VariableSparsity
        if maxdegree === nothing
            error("`maxdegree` cannot be `nothing` when `sparsity` is no `NoSparsity`.")
        end
        c = Certificate.ChordalIdeal(cone, basis, maxdegree)
    elseif sparsity isa MonomialSparsity
        error("Monomial sparsity not implemented yet")
    else
        @assert sparsity isa NoSparsity
        if basis isa AbstractPolynomialBasis
            c = Certificate.FixedBasis(cone, basis)
        else
            c = default_ideal_certificate(domain, cone, basis, newton_polytope, maxdegree)
        end
    end
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
function default_certificate(::BasicSemialgebraicSet, ::VariableSparsity,
                             ideal_certificate::Certificate.ChordalIdeal, cone,
                             basis, maxdegree)
    return Certificate.ChordalPutinar(
         cone, basis, maxdegree)
end
function default_certificate(::BasicSemialgebraicSet, ::MonomialSparsity,
                             ideal_certificate, cone,
                             basis, maxdegree)
    error("Monomial sparsity not implemented yet")
end
function default_certificate(::BasicSemialgebraicSet, ::NoSparsity,
                             ideal_certificate, cone, basis, maxdegree)
    return Certificate.Putinar(
        ideal_certificate, cone, basis, maxdegree)
end
function JuMP.moi_set(
    cone::SOSLikeCone,
    monos::AbstractVector{<:MP.AbstractMonomial};
    domain::AbstractSemialgebraicSet=FullSpace(),
    basis=MonomialBasis,
    newton_polytope::Tuple=tuple(),
    maxdegree::Union{Nothing, Int}=MP.maxdegree(monos),
    sparsity::Sparsity=NoSparsity(),
    newton_of_remainder::Bool=false,
    ideal_certificate=default_ideal_certificate(
        domain, cone, basis, newton_polytope, maxdegree, sparsity, newton_of_remainder),
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
