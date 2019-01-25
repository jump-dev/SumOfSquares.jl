export slack, certificate_monomials, lagrangian_multipliers

function JuMP.moi_set(cone::SOSSubCones,
                      monos::AbstractVector{<:AbstractMonomial};
                      domain::AbstractSemialgebraicSet=FullSpace(),
                      basis=MonomialBasis,
                      newton_polytope::Tuple=tuple(),
                      mindegree=MultivariatePolynomials.mindegree(monos),
                      maxdegree=MultivariatePolynomials.maxdegree(monos))
    return SOSPolynomialSet(domain, cone, basis, monos, tuple(), mindegree,
                            maxdegree)
end

function JuMP.build_constraint(_error::Function, p, cone::SOSSubCones; kws...)
    coefs = PolyJuMP.non_constant_coefficients(p)
    monos = monomials(p)
    set = JuMP.moi_set(cone, monos; kws...)
    shape = PolyJuMP.PolynomialShape(monos)
    return JuMP.VectorConstraint(coefs, set, shape)
end

slack(cref::JuMP.ConstraintRef) =  MOI.get(cref.model, Slack(), cref)

function certificate_monomials(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, CertificateMonomials(), cref)
end

function lagrangian_multipliers(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, LagrangianMultipliers(), cref)
end

function MOI.is_set_by_optimize(::Union{Slack, CertificateMonomials,
                                        LagrangianMultipliers})
    return true
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
                          newton_polytope=(y, newton_polytope...))
end
