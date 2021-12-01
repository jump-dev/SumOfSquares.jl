using DynamicPolynomials
@polyvar x[1:3]
p = sum(x)^2
using SumOfSquares
S = algebraicset([xi^2 - 1 for xi in x])

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

function min_algebraic(S)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = S)
    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
end

min_algebraic(S)

@show S.I.gröbnerbasis
S.I.algo

const MP = MultivariatePolynomials
const SS = SemialgebraicSets
struct HypercubeIdeal{V} <: SS.AbstractPolynomialIdeal
    variables::Vector{V}
end
struct HypercubeSet{V} <: SS.AbstractAlgebraicSet
    ideal::HypercubeIdeal{V}
end
MP.changecoefficienttype(set::HypercubeSet, ::Type) = set
SS.ideal(set::HypercubeSet) = set.ideal
function Base.rem(p, set::HypercubeIdeal)
    return MP.polynomial(map(MP.terms(p)) do term
        mono = MP.monomial(term)
        new_mono = one(mono)
        for (var, exp) in powers(mono)
            if var in set.variables
                exp = rem(exp, 2)
            end
            new_mono *= var^exp
        end
        MP.coefficient(term) * new_mono
    end)
end

H = HypercubeSet(HypercubeIdeal(x))

min_algebraic(H)

function min_algebraic_rational(S, d)
    model = SOSModel(solver)
    @variable(model, q, SOSPoly(MP.monomials(x, 0:d)))
    deno = q + 1
    @constraint(model, c, deno * p >= deno, domain = S)
    optimize!(model)
    @show termination_status(model)
end

min_algebraic_rational(H, 0)

min_algebraic_rational(H, 1)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

