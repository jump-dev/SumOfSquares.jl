export getslack, addpolyeqzeroconstraint, addpolynonnegativeconstraint
import JuMP.getdual, PolyJuMP.getslack

type SOSConstraint
  slack::MatPolynomial{JuMP.Variable}
  lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
  x::MonomialVector
end

function getslack(c::SOSConstraint)
  getvalue(c.slack)
end

function getdual(c::SOSConstraint)
  a = [getdual(lc) for lc in c.lincons]
  PseudoExpectation(a, c.x)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p)
  constraints = [JuMP.constructconstraint!(t.α, :(==)) for t in p]
  JuMP.addVectorizedConstraint(m, constraints)
end

function addpolynonnegativeconstraint(m::JuMP.Model, P::Matrix)
  n = Base.LinAlg.checksquare(P)
  if !issymmetric(P)
    error("The polynomial matrix constrained to be SOS must be symmetric")
  end
  y = polyvecvar(string(gensym()), 1:n)
  p = dot(y, P*y)
  addpolynonnegativeconstraint(m, p)
end
function addpolynonnegativeconstraint(m::JuMP.Model, p)
  Z = getmonomialsforcertificate(p.x)
  slack = createnonnegativepoly(m, :Gram, Z)
  q = p - slack
  lincons = addpolyeqzeroconstraint(m, q)
  SOSConstraint(slack, lincons, q.x)
end
