using JuMP
import JuMP.getvalue
export @SOSvariable, @SOSconstraint, getslack

function getvalue(p::MatPolynomial{JuMP.Variable})
  MatPolynomial(map(getvalue, p.Q), p.x)
end

function addpolyeqzeroconstraint(m::JuMP.Model, p)
  constraints = [JuMP.constructconstraint!(t.α, :(==)) for t in p]
  JuMP.addVectorizedConstraint(m, constraints)
end
function addsosconstraint(m::JuMP.Model, p)
  Z = getmonomialsforcertificate(p.x)
  @show Z
  slack = MatPolynomial{JuMP.Variable}((i,j) -> Variable(m, -Inf, Inf, :Cont), Z)
  push!(m.varCones, (:SDP, slack.Q[1].col:slack.Q[end].col))
  lincons = addpolyeqzeroconstraint(m, p - slack)
  SOSConstraintRef(slack, lincons)
end

macro SOSvariable(args...)
  length(args) <= 1 &&
  variable_error(args, "Expected model as first argument, then variable information.")
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
    variable_error(args, "SOS variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
  elseif isexpr(x,:call)
    explicit_comparison = true
    if x.args[1] == :>= || x.args[1] == :≥
      # x >= lb
      var = x.args[2]
      @assert length(x.args) == 3
      lb = esc_nonconstant(x.args[3])
      if lb != 0
        variable_error(args, "SOS variable declaration does not support the form ... >= lb with nonzero lb.")
      end
      sos = true
    elseif x.args[1] == :<= || x.args[1] == :≤
      # x <= ub
      # NB: May also be lb <= x, which we do not support
      #     We handle this later in the macro
      variable_error(args, "SOS variable declaration does not support the form ... <= ub.")
    elseif x.args[1] == :(==)
      variable_error(args, "SOS variable declaration does not support the form ... == value.")
    else
      # Its a comparsion, but not using <= ... <=
      variable_error(args, "Unexpected syntax $(string(x)).")
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
  @show kwargs
  @show extra

  anonvar = isexpr(var, :vect) || isexpr(var, :vcat)
  anonvar && explicit_comparison && error("Cannot use explicit bounds via >=, <= with an anonymous variable")
  variable = gensym()
  quotvarname = anonvar ? :(:__anon__) : quot(getname(var))
  escvarname  = anonvar ? variable     : esc(getname(var))

  # process keyword arguments
  value = NaN
  obj = nothing
  inconstraints = nothing
  coefficients = nothing
  for ex in kwargs
    kwarg = ex.args[1]
    if kwarg == :start
      value = esc(ex.args[2])
    elseif kwarg == :objective
      obj = esc(ex.args[2])
    elseif kwarg == :inconstraints
      inconstraints = esc(ex.args[2])
    elseif kwarg == :coefficients
      coefficients = esc(ex.args[2])
    elseif kwarg == :basename
      quotvarname = esc(ex.args[2])
    elseif kwarg == :lowerbound
      haslb && variable_error(args, "Cannot specify variable lowerbound twice")
      lb = esc_nonconstant(ex.args[2])
      haslb = true
    elseif kwarg == :upperbound
      hasub && variable_error(args, "Cannot specify variable upperbound twice")
      ub = esc_nonconstant(ex.args[2])
      hasub = true
    else
      variable_error(args, "Unrecognized keyword argument $kwarg")
    end
  end

  if (obj !== nothing || inconstraints !== nothing || coefficients !== nothing) &&
    (obj === nothing || inconstraints === nothing || coefficients === nothing)
    variable_error(args, "Must provide 'objective', 'inconstraints', and 'coefficients' arguments all together for column-wise modeling")
  end

  sdp = any(t -> (t == :SDP), extra)
  symmetric = (sdp || any(t -> (t == :Symmetric), extra))
  extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag

  # Determine variable type (if present).
  # Types: default is continuous (reals)
  if length(extra) > 0
    if t == :Fixed
      variable_error(args, "Unexpected extra arguments when declaring a fixed variable")
    end
    if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
      gottype = true
      t = extra[1]
    end

    if t == :Bin
      if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
        variable_error(args, "Bounds other than [0, 1] may not be specified for binary variables.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
      else
        lb = 0.0
        ub = 1.0
      end
    end

    !gottype && variable_error(args, "Syntax error")
  end

  # Handle the column generation functionality
  if coefficients !== nothing
    !isa(var,Symbol) &&
    variable_error(args, "Can only create one variable at a time when adding to existing constraints.")

    return assert_validmodel(m, quote
      $variable = Variable($m,$lb,$ub,$(quot(t)),$obj,$inconstraints,$coefficients,UTF8String(string($quotvarname)),$value)
      $(anonvar ? variable : :($escvarname = $variable))
    end)
  end

  if isa(var,Symbol)
    # Easy case - a single variable
    sdp && variable_error(args, "Cannot add a semidefinite scalar variable")
    @assert !anonvar
    return assert_validmodel(m, quote
      $variable = Variable($m,$lb,$ub,$(quot(t)),UTF8String(string($quotvarname)),$value)
      registervar($m, $quotvarname, $variable)
      $escvarname = $variable
    end)
  end
  isa(var,Expr) || variable_error(args, "Expected $var to be a variable name")

  # We now build the code to generate the variables (and possibly the JuMPDict
  # to contain them)
  refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var, variable)
  clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? () : idxsets[i])

  code = :( $(refcall) = Variable($m, $lb, $ub, $(quot(t)), EMPTYSTRING, $value) )
  if symmetric
    # Sanity checks on SDP input stuff
    condition == :() ||
    variable_error(args, "Cannot have conditional indexing for SDP variables")
    length(idxvars) == length(idxsets) == 2 ||
    variable_error(args, "SDP variables must be 2-dimensional")
    !symmetric || (length(idxvars) == length(idxsets) == 2) ||
    variable_error(args, "Symmetric variables must be 2-dimensional")
    hasdependentsets(idxvars, idxsets) &&
    variable_error(args, "Cannot have index dependencies in symmetric variables")
    for _rng in idxsets
      isexpr(_rng, :escape) ||
      variable_error(args, "Internal error 1")
      rng = _rng.args[1] # undo escaping
      (isexpr(rng,:(:)) && rng.args[1] == 1 && length(rng.args) == 2) ||
      variable_error(args, "Index sets for SDP variables must be ranges of the form 1:N")
    end

    if !(lb == -Inf && ub == Inf)
      variable_error(args, "Semidefinite or symmetric variables cannot be provided bounds")
    end
    return assert_validmodel(m, quote
      $(esc(idxsets[1].args[1].args[2])) == $(esc(idxsets[2].args[1].args[2])) || error("Cannot construct symmetric variables with nonsquare dimensions")
      (Compat.issymmetric($lb) && Compat.issymmetric($ub)) || error("Bounds on symmetric  variables must be symmetric")
      $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :Variable; lowertri=symmetric))
      $(if sdp
        quote
          push!($(m).varCones, (:SDP, first($variable).col : last($variable).col))
        end
      end)
      push!($(m).dictList, $variable)
      registervar($m, $quotvarname, $variable)
      storecontainerdata($m, $variable, $quotvarname,
      $(Expr(:tuple,idxsets...)),
      $idxpairs, $(quot(condition)))
      $(anonvar ? variable : :($escvarname = $variable))
    end)
  else
    coloncheckcode = Expr(:call,:coloncheck,refcall.args[2:end]...)
    code = :($coloncheckcode; $code)
    return assert_validmodel(m, quote
      $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :Variable))
      isa($variable, JuMPContainer) && pushmeta!($variable, :model, $m)
      push!($(m).dictList, $variable)
      registervar($m, $quotvarname, $variable)
      storecontainerdata($m, $variable, $quotvarname,
      $(Expr(:tuple,map(clear_dependencies,1:length(idxsets))...)),
      $idxpairs, $(quot(condition)))
      $(anonvar ? variable : :($escvarname = $variable))
    end)
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
  JuMP.isexpr(x,:call) && length(x.args) == 3 || error("in @SDconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
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
