export getslack, addpolyeqzeroconstraint, addpolynonnegativeconstraint

type SOSConstraint{C}
    slack::MatPolynomial{C, JuMP.Variable}
    lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
    x::MonomialVector{C}
end

function JuMP.getdual(c::SOSConstraint)
    a = [getdual(lc) for lc in c.lincons]
    Measure(a, c.x)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p, domain::FullSpace)
    constraints = [JuMP.constructconstraint!(AffExpr(t.α), :(==)) for t in p]
    JuMP.addVectorizedConstraint(m, constraints)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p, domain::AlgebraicSet)
    if !isempty(domain.p)
        warn("Equality on algebraic set has not been implemented yet, ignoring the domain")
    end
    addpolyeqzeroconstraint(m, p, FullSpace())
end

function addpolyeqzeroconstraint(m::JuMP.Model, p, domain::BasicSemialgebraicSet)
    addpolynonnegativeconstraint(m,  p, domain)
    addpolynonnegativeconstraint(m, -p, domain)
    nothing
end

function matconstraux{C}(::Type{PolyVar{C}}, m::JuMP.Model, P::Matrix, domain::AbstractBasicSemialgebraicSet)
    n = Base.LinAlg.checksquare(P)
    if !issymmetric(P)
        throw(ArgumentError("The polynomial matrix constrained to be SOS must be symmetric"))
    end
    y = polyvecvar(PolyVar{C}, string(gensym()), 1:n)
    p = dot(y, P*y)
    addpolynonnegativeconstraint(m, p, domain)
end

for T in (FullSpace, AlgebraicSet, BasicSemialgebraicSet)
    @eval begin
        addpolynonnegativeconstraint{T<:VectorOfPolyType{false}}(m::JuMP.Model, P::Matrix{T}, domain::$T) = matconstraux(PolyVar{false}, m, P, domain)
        addpolynonnegativeconstraint{T<:VectorOfPolyType{true}}(m::JuMP.Model, P::Matrix{T}, domain::$T) = matconstraux(PolyVar{true}, m, P, domain)
    end
end

function addpolynonnegativeconstraint(m::JuMP.Model, p, domain::FullSpace)
    # FIXME If p is a MatPolynomial, p.x will not be correct
    Z = getmonomialsforcertificate(p.x)
    slack = createnonnegativepoly(m, Z, :Cont)
    q = p - slack
    lincons = addpolyeqzeroconstraint(m, q, domain)
    SOSConstraint(slack, lincons, q.x)
end

function addpolynonnegativeconstraint(m::JuMP.Model, p, domain::AlgebraicSet)
    if !isempty(domain.p)
        warn("Equality on algebraic set has not been implemented yet, ignoring the domain")
    end
    addpolynonnegativeconstraint(m, p, FullSpace())
end

function addpolynonnegativeconstraint(m::JuMP.Model, p, domain::BasicSemialgebraicSet)
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
    addpolynonnegativeconstraint(m, p, domain.V)
end
