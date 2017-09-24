export DSOSCone, SDSOSCone, SOSCone, getslack, addpolyconstraint!
export CoDSOSCone, CoSDSOSCone, CoSOSCone

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

function addpolyconstraint!(m::JuMP.Model, p, set::NonNegPolySubCones, domain::BasicSemialgebraicSet)
    mindeg, maxdeg = extdegree(p)
    for q in domain.p
        mindegq, maxdegq = extdegree(q)
        mind = mindeg - mindegq
        maxd = maxdeg - maxdegq
        mind = max(0, Int(floor(mind / 2)))
        maxd = Int(ceil(maxd / 2))
        # FIXME handle the case where `p`, `q_i`, ...  do not have the same variables
        # so instead of `variable(p)` we would have the union of them all
        @assert variables(q) ⊆ variables(p)
        s = createpoly(m, _varconetype(set)(monomials(variables(p), mind:maxd)), :Cont)
        p -= s*q
    end
    addpolyconstraint!(m, p, set, domain.V)
end
