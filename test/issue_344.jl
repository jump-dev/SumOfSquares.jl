using JuMP
using SumOfSquares
using DynamicPolynomials
function sos_min(sparsity, d, obj, dom)
    model =
        Model(() -> MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}()))
    @variable(model, t)
    @objective(model, Min, t)
    con_ref = @constraint(
        model,
        obj - t in SOSCone(),
        sparsity = sparsity,
        maxdegree = d,
        domain = dom
    )
    return optimize!(model)
    #return moment_matrix(con_ref)
end
@polyvar x[1:3] l[1:3]
my_game_obj = x[1]^2 + x[2]^2 + 1.0 * x[3]^2
my_eq_list = [
    -1.0 + l[3] - 3l[2] - 3x[2]^2 - x[1] * x[2],
    3x[2] * l[2] - x[2]^2 * l[2] - x[1]^2 * l[2],
    l[3] - x[3] * l[3],
    -2 + x[3] + 2x[1] + 2x[1] * l[1],
    -2 - x[3] + 2x[2] + 2x[2] * l[1],
    2l[1] - x[3] * l[1] - x[2]^2 * l[1] - x[1]^2 * l[1],
]
my_ineq_list = [
    3x[3] - x[2]^2 - x[1]^2,
    1 - x[3],
    l[2],
    l[3],
    2 - x[3] - x[2]^2 - x[1]^2,
    l[1],
]
my_Dom = basicsemialgebraicset(algebraicset(my_eq_list), my_ineq_list)
ν = sos_min(Sparsity.Monomial(ChordalCompletion()), 4, my_game_obj, my_Dom)
