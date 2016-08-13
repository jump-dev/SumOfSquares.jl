using JuMP
export @SOSconstraint

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


function addpolyeqzeroconstraint(m::JuMP.Model, p)
  constraints = [JuMP.constructconstraint!(t.α, :(==)) for t in p]
  JuMP.addVectorizedConstraint(m, constraints)
end
function addsosconstraint(m::JuMP.Model, p)
  Z = getmonomialsforcertificate(p.x)
  slack = MatPolynomial{JuMP.Variable}((i,j) -> Variable(m, -Inf, Inf, :Cont), Z)
  push!(m.varCones, (:SDP, slack.Q[1].col:slack.Q[end].col))
  addpolyeqzeroconstraint(m, p - slack)
end
