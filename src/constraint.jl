export getslack, addpolyeqzeroconstraint, addpolynonnegativeconstraint
import JuMP.getdual, PolyJuMP.getslack

type SOSConstraint{C}
    slack::MatPolynomial{C, JuMP.Variable}
    lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
    x::MonomialVector{C}
end

function getslack(c::SOSConstraint)
    getvalue(c.slack)
end

function getdual(c::SOSConstraint)
    a = [getdual(lc) for lc in c.lincons]
    Measure(a, c.x)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p, domain::AlgebraicSet)
    @assert isempty(domain.p)
    constraints = [JuMP.constructconstraint!(t.α, :(==)) for t in p]
    JuMP.addVectorizedConstraint(m, constraints)
end

function matconstraux{C}(::Type{PolyVar{C}}, m::JuMP.Model, P::Matrix, domain::BasicSemialgebraicSet)
    n = Base.LinAlg.checksquare(P)
    if !issymmetric(P)
        throw(ArgumentError("The polynomial matrix constrained to be SOS must be symmetric"))
    end
    y = polyvecvar(PolyVar{C}, string(gensym()), 1:n)
    p = dot(y, P*y)
    addpolynonnegativeconstraint(m, p, domain)
end
addpolynonnegativeconstraint{T<:VectorOfPolyType{false}}(m::JuMP.Model, P::Matrix{T}, domain::BasicSemialgebraicSet) = matconstraux(PolyVar{false}, m, P, domain)
addpolynonnegativeconstraint{T<:VectorOfPolyType{true}}(m::JuMP.Model, P::Matrix{T}, domain::BasicSemialgebraicSet) = matconstraux(PolyVar{true}, m, P, domain)

function addpolynonnegativeconstraint(m::JuMP.Model, p, domain::BasicSemialgebraicSet)
    # TODO We might want to add this as a function in MultivariatePolynomials.jl
    mindeg, maxdeg = extdeg(p)
    for q in domain.p
        mindegq, maxdegq = extdeg(q)
        mind = mindeg - mindegq
        maxd = maxdeg - maxdegq
        mind = max(0, Int(floor(mind / 2)))
        maxd = Int(ceil(maxd / 2))
        # FIXME handle the case where `p`, `q_i`, ...  do not have the same variables
        # so instead of `var(p)` we would have the union of them all
        @assert vars(q) == vars(p)
        s = createnonnegativepoly(m, :Gram,  MonomialVector(vars(p), mind:maxd))
        p -= s*q
    end
    Z = getmonomialsforcertificate(p.x)
    slack = createnonnegativepoly(m, :Gram, Z)
    q = p - slack
    lincons = addpolyeqzeroconstraint(m, q, domain.V)
    SOSConstraint(slack, lincons, q.x)
end
