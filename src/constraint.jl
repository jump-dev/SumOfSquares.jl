export DSOSCone, SDSOSCone, SOSCone
export CoDSOSCone, CoSDSOSCone, CoSOSCone
export getslack, certificate_monomials, addpolyconstraint!

struct DSOSCone end
struct CoDSOSCone end
_varconetype(::DSOSCone) = DSOSPoly
_nococone(::CoDSOSCone) = DSOSCone()

struct SDSOSCone end
struct CoSDSOSCone end
_varconetype(::SDSOSCone) = SDSOSPoly
_nococone(::CoSDSOSCone) = SDSOSCone()

struct SOSCone end
struct CoSOSCone end
_varconetype(::SOSCone) = SOSPoly
_nococone(::CoSOSCone) = SOSCone()

_varconetype(::NonNegPoly) = Poly{true}

const SOSLikeCones = Union{DSOSCone, SDSOSCone, SOSCone, NonNegPoly}
const CoSOSLikeCones = Union{CoDSOSCone, CoSDSOSCone, CoSOSCone}
const NonNegPolySubCones = Union{CoSOSLikeCones, SOSLikeCones}

struct SOSConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JS<:JuMP.AbstractJuMPScalar, JC<:JuMP.AbstractConstraint}
    # JS is AffExpr for CoSOS and is Variable for SOS
    slack::MatPolynomial{JS, MT, MVT}
    lincons::Vector{JuMP.ConstraintRef{JuMP.Model, JC}}
    x::MVT
end

certificate_monomials(c::SOSConstraint) = c.slack.x

function JuMP.getdual(c::SOSConstraint)
    a = [getdual(lc) for lc in c.lincons]
    Measure(a, c.x)
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace)
    constraints = JuMP.constructconstraint!.(coefficients(p), :(==))
    JuMP.addVectorizedConstraint(m, constraints)
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::AbstractAlgebraicSet)
    addpolyconstraint!(m, rem(p, ideal(domain)), s, FullSpace())
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::BasicSemialgebraicSet)
    addpolyconstraint!(m,  p, NonNegPoly(), domain)
    addpolyconstraint!(m, -p, NonNegPoly(), domain)
    nothing
end

function addpolyconstraint!(m::JuMP.Model, P::Matrix{PT}, ::PSDCone, domain::AbstractBasicSemialgebraicSet) where PT <: APL
    n = Base.LinAlg.checksquare(P)
    if !issymmetric(P)
        throw(ArgumentError("The polynomial matrix constrained to be SOS must be symmetric"))
    end
    y = [similarvariable(PT, gensym()) for i in 1:n]
    p = dot(y, P * y)
    addpolyconstraint!(m, p, NonNegPoly(), domain)
end

function _createslack(m, x, set::SOSLikeCones)
    createpoly(m, _varconetype(set)(x), :Cont)
end
function _matposynomial(m, x)
    p = _matpolynomial(m, x, :Cont)
    m.colLower[map(q -> q.col, p.Q)] = 0.
    p
end
function _createslack(m, x, set::CoSOSLikeCones)
    _matplus(_createslack(m, x, _nococone(set)), _matposynomial(m, x))
end

function addpolyconstraint!(m::JuMP.Model, p, set::NonNegPolySubCones, domain::AbstractAlgebraicSet)
    r = rem(p, ideal(domain))
    X = getmonomialsforcertificate(monomials(r))
    slack = _createslack(m, X, set)
    q = r - slack
    lincons = addpolyconstraint!(m, q, ZeroPoly(), domain)
    SOSConstraint(slack, lincons, monomials(q))
end

function addpolyconstraint!(m::JuMP.Model, p, set::NonNegPolySubCones, domain::BasicSemialgebraicSet;
                            mindegree=MultivariatePolynomials.mindegree(p),
                            maxdegree=MultivariatePolynomials.maxdegree(p))
    for q in domain.p
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
        s2 = createpoly(m, _varconetype(set)(monomials(variables(p), mindegree_s:maxdegree_s)), :Cont)
        p -= s2 * q
    end
    addpolyconstraint!(m, p, set, domain.V)
end
