import TypedPolynomials
TypedPolynomials.@polyvar x y
using SumOfSquares
import CSDP
model = SOSModel(CSDP.Optimizer)
con_ref = @constraint(model, 2x^4 + 5y^4 - x^2*y^2 >= -2(x^3*y + x + 1))
optimize!(model)
solution_summary(model)

sos_decomposition(con_ref)

using BenchmarkTools
@btime let
    TypedPolynomials.@polyvar x y
    S = @set x^3 * y + x == 2x^2 * y^2 && 3x^4 == y
    SemialgebraicSets.compute_gröbner_basis!(S.I)
end

import DynamicPolynomials
@btime let
    DynamicPolynomials.@polyvar x y
    S = @set x^3 * y + x == 2x^2 * y^2 && 3x^4 == y
    SemialgebraicSets.compute_gröbner_basis!(S.I)
end

TypedPolynomials.@polyvar x y
S = @set x^3 * y + x == 2x^2 * y^2 && 3x^4 == y
poly = -6x - 4y^3 + 2x*y^2 + 6x^3 - 3y^4 + 13x^2 * y^2
model = Model(CSDP.Optimizer)
con_ref = @constraint(model, poly in SOSCone(), domain = S)
optimize!(model)
solution_summary(model)

dec = sos_decomposition(con_ref, 1e-6)

rem(dec - poly, S.I)

rem(gram_matrix(con_ref) - poly, S.I)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
