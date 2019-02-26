export certificate_monomials, gram_matrix, lagrangian_multipliers
export SOSMatrixCone, SOSConvexCone

function JuMP.reshape_set(set::SOSPolynomialSet, ::PolyJuMP.PolynomialShape)
    return set.cone
end
function JuMP.moi_set(cone::SOSLikeCone,
                      monos::AbstractVector{<:AbstractMonomial};
                      domain::AbstractSemialgebraicSet=FullSpace(),
                      basis=MonomialBasis,
                      newton_polytope::Tuple=tuple(),
                      mindegree=MultivariatePolynomials.mindegree(monos),
                      maxdegree=MultivariatePolynomials.maxdegree(monos))
    return SOSPolynomialSet(domain, cone, basis, monos, newton_polytope,
                            mindegree, maxdegree)
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:SOSPolynomialSet{<:AbstractAlgebraicSet}})
    return [SOSPolynomialBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:DiagonallyDominantConeTriangle})
    return [DiagonallyDominantBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:ScaledDiagonallyDominantConeTriangle})
    return [ScaledDiagonallyDominantBridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{PositiveSemidefinite2x2ConeTriangle})
    return [PositiveSemidefinite2x2Bridge]
end

function PolyJuMP.bridges(::Type{<:MOI.AbstractVectorFunction},
                          ::Type{<:SOSPolynomialSet{<:BasicSemialgebraicSet}})
    return [SOSPolynomialInSemialgebraicSetBridge]
end

function JuMP.build_constraint(_error::Function, p, cone::SOSLikeCone; kws...)
    coefs = PolyJuMP.non_constant_coefficients(p)
    monos = monomials(p)
    set = JuMP.moi_set(cone, monos; kws...)
    shape = PolyJuMP.PolynomialShape(monos)
    return PolyJuMP.bridgeable(JuMP.VectorConstraint(coefs, set, shape),
                               JuMP.moi_function_type(typeof(coefs)),
                               typeof(set))
end

function gram_matrix(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, GramMatrixAttribute(), cref)
end

function MultivariateMoments.moment_matrix(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, MomentMatrixAttribute(), cref)
end

# Equivalent but more efficient than moment_matrix(cref).x as it does not need
# to query any result from the solver
function certificate_monomials(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, CertificateMonomials(), cref)
end

function lagrangian_multipliers(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, LagrangianMultipliers(), cref)
end

struct SOSMatrixCone <: PolyJuMP.PolynomialSet end

function JuMP.build_constraint(_error::Function, P::Matrix{PT},
                               ::SOSMatrixCone;
                               newton_polytope::Tuple = tuple(),
                               kws...) where PT <: APL
    n = LinearAlgebra.checksquare(P)
    if !issymmetric(P)
        _error("The polynomial matrix constrained to be SOS must be symmetric.")
    end
    y = [similarvariable(PT, gensym()) for i in 1:n]
    p = dot(y, P * y)
    # TODO Newton_polytope=(y,) may not be the best idea if exact newton
    #      polytope computation is used.
    #      See "Sum-of-Squares Matrices" notebook
    JuMP.build_constraint(_error, p, SOSCone();
                          newton_polytope=(y, newton_polytope...), kws...)
end

struct SOSConvexCone <: PolyJuMP.PolynomialSet end

function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               ::SOSConvexCone; kws...)
    hessian = differentiate(p, variables(p), 2)
    JuMP.build_constraint(_error, hessian, SOSMatrixCone(); kws...)
end
