export getslack, addpolyeqzeroconstraint, addpolynonnegativeconstraint
import JuMP.getdual

type SOSConstraintRef
  slack::MatPolynomial{JuMP.Variable}
  lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
  x::MonomialVector
end

function getslack(c::SOSConstraintRef)
  getvalue(c.slack)
end

function getdual(c::SOSConstraintRef)
  a = [getdual(lc) for lc in c.lincons]
  PseudoExpectation(a, c.x)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p)
  constraints = [JuMP.constructconstraint!(t.α, :(==)) for t in p]
  JuMP.addVectorizedConstraint(m, constraints)
end

function addpolynonnegativeconstraint(m::JuMP.Model, sos::SOS, P::Matrix)
  n = Base.LinAlg.checksquare(P)
  if !issymmetric(P)
    error("The polynomial matrix constrained to be SOS must be symmetric")
  end
  y = polyvecvar(string(gensym()), 1:n)
  p = dot(y, P*y)
  addpolynonnegativeconstraint(m, sos, p)
end
function addpolynonnegativeconstraint(m::JuMP.Model, sos::SOS, p)
  Z = getmonomialsforcertificate(p.x)
  slack = createnonnegativepoly(m, sos, :Gram, Z)
  q = p - slack
  lincons = addpolyeqzeroconstraint(m, q)
  SOSConstraintRef(slack, lincons, q.x)
end
