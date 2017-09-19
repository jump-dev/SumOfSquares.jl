export SOSCone, getslack, addpolyconstraint!

struct SOSCone end
const NonNegPolySubCones = Union{NonNegPoly, SOSCone}

struct SOSConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}}
    slack::MatPolynomial{JuMP.Variable, MT, MVT}
    lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
    x::MVT
end

function JuMP.getdual(c::SOSConstraint)
    a = [getdual(lc) for lc in c.lincons]
    Measure(a, c.x)
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace)
    constraints = JuMP.constructconstraint!.(AffExpr.(coefficients(p)), :(==))
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

function addpolyconstraint!(m::JuMP.Model, p, ::NonNegPolySubCones, domain::AbstractAlgebraicSet)
    r = rem(p, ideal(domain))
    X = getmonomialsforcertificate(monomials(r))
    slack = createpoly(m, Poly{true}(X), :Cont)
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
        s = createpoly(m, Poly{true}(monomials(variables(p), mind:maxd)), :Cont)
        p -= s*q
    end
    addpolyconstraint!(m, p, set, domain.V)
end
