export SOSCone, getslack, addpolyconstraint!

type SOSCone end
const NonNegPolySubCones = Union{NonNegPoly, SOSCone}

type SOSConstraint{C}
    slack::MatPolynomial{C, JuMP.Variable}
    lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
    x::MonomialVector{C}
end

function JuMP.getdual(c::SOSConstraint)
    a = [getdual(lc) for lc in c.lincons]
    Measure(a, c.x)
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace)
    constraints = [JuMP.constructconstraint!(AffExpr(t.α), :(==)) for t in p]
    JuMP.addVectorizedConstraint(m, constraints)
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::AlgebraicSet)
    if !isempty(domain.p)
        warn("Equality on algebraic set has not been implemented yet, ignoring the domain")
    end
    addpolyconstraint!(m, p, FullSpace())
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::BasicSemialgebraicSet)
    addpolyconstraint!(m,  p, NonNegPoly(), domain)
    addpolyconstraint!(m, -p, NonNegPoly(), domain)
    nothing
end

function matconstraux{C}(::Type{PolyVar{C}}, m::JuMP.Model, P::Matrix, domain::AbstractBasicSemialgebraicSet)
    n = Base.LinAlg.checksquare(P)
    if !issymmetric(P)
        throw(ArgumentError("The polynomial matrix constrained to be SOS must be symmetric"))
    end
    y = polyvecvar(PolyVar{C}, string(gensym()), 1:n)
    p = dot(y, P*y)
    addpolyconstraint!(m, p, NonNegPoly(), domain)
end

for T in (FullSpace, AlgebraicSet, BasicSemialgebraicSet)
    @eval begin
        addpolyconstraint!{T<:VectorOfPolyType{false}}(m::JuMP.Model, P::Matrix{T}, ::PSDCone, domain::$T) = matconstraux(PolyVar{false}, m, P, domain)
        addpolyconstraint!{T<:VectorOfPolyType{true}}(m::JuMP.Model, P::Matrix{T}, ::PSDCone, domain::$T) = matconstraux(PolyVar{true}, m, P, domain)
    end
end

function addpolyconstraint!(m::JuMP.Model, p, ::Union{NonNegPoly, SOSCone}, domain::FullSpace)
    # FIXME If p is a MatPolynomial, p.x will not be correct
    Z = getmonomialsforcertificate(p.x)
    slack = createnonnegativepoly(m, Z, :Cont)
    q = p - slack
    lincons = addpolyconstraint!(m, q, ZeroPoly(), domain)
    SOSConstraint(slack, lincons, q.x)
end

function addpolyconstraint!(m::JuMP.Model, p, s::NonNegPolySubCones, domain::AlgebraicSet)
    if !isempty(domain.p)
        warn("Equality on algebraic set has not been implemented yet, ignoring the domain")
    end
    addpolyconstraint!(m, p, s, FullSpace())
end

function addpolyconstraint!(m::JuMP.Model, p, set::NonNegPolySubCones, domain::BasicSemialgebraicSet)
    mindeg, maxdeg = extdeg(p)
    for q in domain.p
        mindegq, maxdegq = extdeg(q)
        mind = mindeg - mindegq
        maxd = maxdeg - maxdegq
        mind = max(0, Int(floor(mind / 2)))
        maxd = Int(ceil(maxd / 2))
        # FIXME handle the case where `p`, `q_i`, ...  do not have the same variables
        # so instead of `var(p)` we would have the union of them all
        @assert vars(q) ⊆ vars(p)
        s = createnonnegativepoly(m, MonomialVector(vars(p), mind:maxd), :Cont)
        p -= s*q
    end
    addpolyconstraint!(m, p, set, domain.V)
end
