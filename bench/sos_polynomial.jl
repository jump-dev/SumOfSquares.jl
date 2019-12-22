using SumOfSquares
using DynamicPolynomials

function sos_polynomial(n)
    @polyvar x
    monos = monomials(x, 0:2n)
    set = SumOfSquares.SOSPolynomialSet(
        SumOfSquares.FullSpace(), monos,
        SumOfSquares.Certificate.Remainder(
            SumOfSquares.SOSCone(),
            SumOfSquares.MonomialBasis,
            tuple()
        )
    )
    model = MOIU.UniversalFallback(MOIU.Model{Float64}())
    x = MOI.add_variables(model, n)
    scalar = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0)
    func = MOI.Utilities.vectorize([scalar for i in eachindex(monos)])
    BT = MOI.Bridges.Constraint.concrete_bridge_type(
        SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Float64},
        typeof(func),
        typeof(set))
    @time MOI.Bridges.Constraint.bridge_constraint(BT, model, func, set)
    return
end
#sos_polynomial(200)
sos_polynomial(2)
