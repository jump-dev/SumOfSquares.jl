using JuMP
import JuMP: getvalue, validmodel, addtoexpr_reorder
using Base.Meta

export @SOSvariable, @SOSconstraint, getslack

function freshmatpoly(m::JuMP.Model, Z::Union{MonomialVector,Vector})
  MatPolynomial{JuMP.Variable}((i,j) -> Variable(m, -Inf, Inf, :Cont), Z)
end

function freshvecpoly(m::JuMP.Model, Z::Union{MonomialVector,Vector})
  VecPolynomial{JuMP.Variable}((i) -> Variable(m, -Inf, Inf, :Cont), Z)
end

function getvalue(p::VecPolynomial{JuMP.Variable})
  VecPolynomial(map(getvalue, p.a), p.x)
end
function getvalue(p::MatPolynomial{JuMP.Variable})
  MatPolynomial(map(getvalue, p.Q), p.x)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p)
  constraints = [JuMP.constructconstraint!(t.α, :(==)) for t in p]
  JuMP.addVectorizedConstraint(m, constraints)
end
function addsosconstraint(m::JuMP.Model, p)
  Z = getmonomialsforcertificate(p.x)
  slack = freshmatpoly(m, Z)
  push!(m.varCones, (:SDP, slack.Q[1].col:slack.Q[end].col))
  lincons = addpolyeqzeroconstraint(m, p - slack)
  SOSConstraintRef(slack, lincons)
end

macro SOSvariable(args...)
  length(args) <= 1 &&
  JuMP.variable_error(args, "Expected model as first argument, then variable information.")
  m = esc(args[1])
  x = args[2]
  extra = vcat(args[3:end]...)

  t = :Cont
  gottype = false
  haslb = false
  hasub = false
  # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
  # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
  if VERSION < v"0.5.0-dev+3231"
    x = comparison_to_call(x)
  end
  explicit_comparison = false
  if isexpr(x,:comparison) # two-sided
    JuMP.variable_error(args, "SOS variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
  elseif isexpr(x,:call)
    explicit_comparison = true
    if x.args[1] == :>= || x.args[1] == :≥
      # x >= lb
      var = x.args[2]
      @assert length(x.args) == 3
      lb = JuMP.esc_nonconstant(x.args[3])
      if lb != 0
        JuMP.variable_error(args, "SOS variable declaration does not support the form ... >= lb with nonzero lb.")
      end
      sos = true
    elseif x.args[1] == :<= || x.args[1] == :≤
      # x <= ub
      # NB: May also be lb <= x, which we do not support
      #     We handle this later in the macro
      JuMP.variable_error(args, "SOS variable declaration does not support the form ... <= ub.")
    elseif x.args[1] == :(==)
      JuMP.variable_error(args, "SOS variable declaration does not support the form ... == value.")
    else
      # Its a comparsion, but not using <= ... <=
      JuMP.variable_error(args, "Unexpected syntax $(string(x)).")
    end
  else
    # No bounds provided - free variable
    # If it isn't, e.g. something odd like f(x), we'll handle later
    var = x
    sos = false
  end

  # separate out keyword arguments
  kwargs = filter(ex->isexpr(ex,:kw), extra)
  extra = filter(ex->!isexpr(ex,:kw), extra)

  variable = gensym()
  quotvarname = quot(getname(var))
  escvarname  = esc(getname(var))

  monotype = :None

  # process keyword arguments
  gram = false
  for ex in kwargs
    kwarg = ex.args[1]
    if kwarg == :grammonomials
      if monotype != :None
        error("Monomials given twice")
      end
      monotype = :Gram
      x = esc(ex.args[2])
    elseif kwarg == :monomials
      if monotype != :None
        error("Monomials given twice")
      end
      monotype = :Classic
      x = esc(ex.args[2])
    else
      JuMP.variable_error(args, "Unrecognized keyword argument $kwarg")
    end
  end

  # Determine variable type (if present).
  # Types: default is continuous (reals)
  if isempty(extra)
    JuMP.variable_error(args, "Missing monomial vector")
  elseif length(extra) > 1
    JuMP.variable_error(args, "Too many extra argument: only expected monomial vector")
  else
    if monotype != :None
      error("Monomials given twice")
    end
    monotype = sos ? :Gram : :Classic
    x = esc(extra[1])
  end
  Z = gensym()

  if monotype == :None
    error("Monomials not given")
  end

  gram = monotype == :Gram

  if isa(var,Symbol)
    # Easy case - a single variable
    return JuMP.assert_validmodel(m, quote
      if $sos
        if $gram
          $Z = $x
        else
          $Z = getmonomialsforcertificate($x)
        end
        # The coefficients of a monomial not in Z do not all have to be zero, only their sum
        $variable = freshmatpoly($m, $Z)
        if !$gram
          addpolyeqzeroconstraint($m, removemonomials(VecPolynomial($variable), $x))
        end
      else
        if $gram
          # We do not want to enforce p(x) = Z(x)^T Q Z(x) to be SOS
          # So we can just compute the monomials of p(x)
          $Z = dot($x, $x).x
        else
          $Z = $x
        end
        $variable = freshvecpoly($m, $Z)
      end
      $escvarname = $variable
    end)
  else
    JuMP.variable_error(args, "Invalid syntax for variable name: $(string(var))")
  end
end

type SOSConstraintRef
  slack::MatPolynomial{JuMP.Variable}
  lincons::Vector{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}}}
end

function getslack(c::SOSConstraintRef)
  getvalue(c.slack)
end

macro SOSconstraint(m, x)
  m = esc(m)

  if isa(x, Symbol)
    error("in @SDConstraint: Incomplete constraint specification $x. Are you missing a comparison (<= or >=)?")
  end

  (x.head == :block) &&
  error("Code block passed as constraint.")
  if VERSION < v"0.5.0-dev+3231"
    x = comparison_to_call(x)
  end
  isexpr(x,:call) && length(x.args) == 3 || error("in @SDconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
  "       expr1 <= expr2\n" * "       expr1 >= expr2")
  # Build the constraint
  # Simple comparison - move everything to the LHS
  sense = x.args[1]
  if sense == :⪰
    sense = :(>=)
  elseif sense == :⪯
    sense = :(<=)
  end
  sense,_ = JuMP._canonicalize_sense(sense)
  lhs = :()
  if sense == :(>=) || sense == :(==)
    lhs = :($(x.args[2]) - $(x.args[3]))
  elseif sense == :(<=)
    lhs = :($(x.args[3]) - $(x.args[2]))
  else
    error("Invalid sense $sense in SOS constraint")
  end
  newaff, parsecode = JuMP.parseExprToplevel(lhs, :q)
  JuMP.assert_validmodel(m, quote
    q = zero(AffExpr)
    $parsecode
    if $sense == :(==)
      addpolyeqzeroconstraint($m, $newaff)
    else
      addsosconstraint($m, $newaff)
    end
  end)
end
