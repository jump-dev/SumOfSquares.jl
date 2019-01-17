function JuMP.build_constraint(_error::Function, p, cone::SOSSubCones;
                               domain::AbstractSemialgebraicSet=FullSpace(),
                               basis=MonomialBasis,
                               newton_polytope::Tuple=tuple(),
                               mindegree=MultivariatePolynomials.mindegree(p),
                               maxdegree=MultivariatePolynomials.maxdegree(p))
    set = SOSPolynomialSet(domain, cone, basis, monomials(p), newton_polytope,
                           mindegree, maxdegree)
    return JuMP.build_constraint(_error, coefficients(p), set)
end

struct CertificateMonomials <: MOI.AbstractConstraintAttribute end
function certificate_monomials(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, CertificateMonomials(), cref)
end

struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute end
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
        _error("The polynomial matrix constrained to be SOS must be symmetric")
    end
    y = [similarvariable(PT, gensym()) for i in 1:n]
    p = dot(y, P * y)
    # TODO Newton_polytope=(y,) may not be the best idea if exact newton
    #      polytope computation is used.
    #      See "Sum-of-Squares Matrices" notebook
    JuMP.build_constraint(_error, p, SOSCone();
                          newton_polytope=(y, newton_polytope...))
end
