using MathOptInterface
const MOI = MathOptInterface

export DSOSCone, SDSOSCone, SOSCone

export CoDSOSCone, CoSDSOSCone, CoSOSCone
export SOSMatrixCone
export getslack, certificate_monomials, lagrangian_multipliers

struct DSOSCone <: PolyJuMP.PolynomialSet end
struct CoDSOSCone <: PolyJuMP.PolynomialSet end
_varconetype(::DSOSCone) = DSOSPoly
_nococone(::CoDSOSCone) = DSOSCone()

struct SDSOSCone <: PolyJuMP.PolynomialSet end
struct CoSDSOSCone <: PolyJuMP.PolynomialSet end
_varconetype(::SDSOSCone) = SDSOSPoly
_nococone(::CoSDSOSCone) = SDSOSCone()

struct SOSCone <: PolyJuMP.PolynomialSet end
struct CoSOSCone <: PolyJuMP.PolynomialSet end
_varconetype(::SOSCone) = SOSPoly
_nococone(::CoSOSCone) = SOSCone()

struct SOSMatrixCone <: PolyJuMP.PolynomialSet end

const SOSLikeCones = Union{DSOSCone, SDSOSCone, SOSCone}
const CoSOSLikeCones = Union{CoDSOSCone, CoSDSOSCone, CoSOSCone}
const SOSSubCones = Union{CoSOSLikeCones, SOSLikeCones}

mutable struct SOSConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JS<:JuMP.AbstractJuMPScalar, F<:MOI.AbstractVectorFunction} <: PolyJuMP.ConstraintDelegate
    # JS is AffExpr for CoSOS and is Variable for SOS
    slack::MatPolynomial{JS, MT, MVT}
    zero_constraint::PolyJuMP.ZeroConstraint{MT, MVT, F}
    lagrangian_multipliers::Vector{MatPolynomial{JS, MT, MVT}}
end
function SOSConstraint(slack::MatPolynomial{JS, MT, MVT},
                       zero_constraint::PolyJuMP.ZeroConstraint{MT, MVT, F}) where {MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JS<:JuMP.AbstractJuMPScalar, F<:MOI.AbstractVectorFunction}
    return SOSConstraint(slack, zero_constraint, MatPolynomial{JS, MT, MVT}[])
end

certificate_monomials(c::PolyJuMP.PolyConstraintRef) = certificate_monomials(PolyJuMP.getdelegate(c))
certificate_monomials(c::SOSConstraint) = c.slack.x
lagrangian_multipliers(c::PolyJuMP.PolyConstraintRef) = lagrangian_multipliers(PolyJuMP.getdelegate(c))
lagrangian_multipliers(c::SOSConstraint) = c.lagrangian_multipliers

JuMP.dual(c::SOSConstraint) = JuMP.dual(c.zero_constraint)

PolyJuMP.getslack(c::SOSConstraint) = JuMP.value(c.slack)

function PolyJuMP.addpolyconstraint!(m::JuMP.Model, P::Matrix{PT}, ::SOSMatrixCone, domain::AbstractBasicSemialgebraicSet, basis) where PT <: APL
    n = LinearAlgebra.checksquare(P)
    if !issymmetric(P)
        throw(ArgumentError("The polynomial matrix constrained to be SOS must be symmetric"))
    end
    y = [similarvariable(PT, gensym()) for i in 1:n]
    p = dot(y, P * y)
    # TODO Newton_polytope=(y,) may not be the best idea if exact newton polytope computation is used.
    #      See "Sum-of-Squares Matrices" notebook
    PolyJuMP.addpolyconstraint!(m, p, SOSCone(), domain, basis; newton_polytope=(y,))
end

function _createslack(model, x, set::SOSLikeCones)
    JuMP.add_variable(model,
                      PolyJuMP.Variable(_varconetype(set)(x), false, false))
end
function _matposynomial(m, x)
    p = _matpolynomial(m, x, false, false)
    for q in p.Q
        JuMP.set_lower_bound(q, 0)
    end
    p
end
function _createslack(m, x, set::CoSOSLikeCones)
    _matplus(_createslack(m, x, _nococone(set)), _matposynomial(m, x))
end

function PolyJuMP.addpolyconstraint!(m::JuMP.Model, p, set::SOSSubCones, domain::AbstractAlgebraicSet, basis; newton_polytope=tuple())
    r = rem(p, ideal(domain))
    X = monomials_half_newton_polytope(monomials(r), newton_polytope)
    slack = _createslack(m, X, set)
    q = r - slack
    zero_constraint = PolyJuMP.addpolyconstraint!(m, q, ZeroPoly(), domain, basis)
    SOSConstraint(slack, zero_constraint)
end

function lagrangian_multiplier(model::JuMP.Model, p, set::SOSSubCones, q, mindegree::Integer, maxdegree::Integer)
    mindegree_q, maxdegree_q = extdegree(q)
    # extdegree's that s^2 should have so that s^2 * p has degrees between mindegree and maxdegree
    mindegree_s2 = mindegree - mindegree_q
    maxdegree_s2 = maxdegree - maxdegree_q
    # extdegree's for s
    mindegree_s = max(0, div(mindegree_s2, 2))
    # If maxdegree_s2 is odd, div(maxdegree_s2,2) would make s^2 have degree up to maxdegree_s2-1
    # for this reason, we take div(maxdegree_s2+1,2) so that s^2 have degree up to maxdegree_s2+1
    maxdegree_s = div(maxdegree_s2 + 1, 2)
    # FIXME handle the case where `p`, `q_i`, ...  do not have the same variables
    # so instead of `variable(p)` we would have the union of them all
    @assert variables(q) ⊆ variables(p)
    monos = monomials(variables(p), mindegree_s:maxdegree_s)
    return JuMP.add_variable(model, PolyJuMP.Variable(_varconetype(set)(monos),
                                                      false, false))
end

function PolyJuMP.addpolyconstraint!(m::JuMP.Model, p, set::SOSSubCones, domain::BasicSemialgebraicSet, basis;
                                     mindegree=MultivariatePolynomials.mindegree(p),
                                     maxdegree=MultivariatePolynomials.maxdegree(p), kws...)
    λ = lagrangian_multiplier.(Ref(m), p, Ref(set), domain.p, mindegree, maxdegree)
    p -= dot(λ, domain.p)
    constraint = PolyJuMP.addpolyconstraint!(m, p, set, domain.V, basis; kws...)
    constraint.lagrangian_multipliers = λ
    return constraint
end
