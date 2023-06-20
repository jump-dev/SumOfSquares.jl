# # Hypercube

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extension/hypercube.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extension/hypercube.ipynb)
# **Contributed by**: Benoît Legat

# Given a Sum-of-Squares constraint on an algebraic set:
# ```math
# g_1(x) = , \ldots, g_m(x) = 0 \Rightarrow p(x) \ge 0.
# ```
# We can either use the certificate:
# ```math
# p(x) = s(x) + \lambda_1(x) g_1(x) + \cdots + \lambda_m(x) g_m(x), s_0(x) \text{ is SOS},
# ```
# or
# ```math
# p(x) \equiv s(x) \pmod{\la g_1(x), \ldots, g_m(x) \ra}, s_0(x) \text{ is SOS}.
# ```
# the second one leads to a *simpler* SDP but needs to compute a *Gr\"obner* basis:
#  * SemialgebraicSets implements Buchberger's algorithm.
#  * The `@set` macro recognizes variable fix, e.g., `x = 1`
#     and provides shortcut.
#  * If you know a \alert{better} way to take modulo,
#    better create your \alert{own} type of algebraic set!

# We illustrate this in this example.

using DynamicPolynomials
@polyvar x[1:3]
p = sum(x)^2
using SumOfSquares
S = algebraicset([xi^2 - 1 for xi in x])

# We will now search for the minimum of `x` over `S` using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.8/installation/#Supported-solvers) for a list of the available choices.

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

# Note that the minimum is in fact `1`.
# Indeed, since each variables is odd (it is either `-1` or `1`)
# and there is an odd number of variables, their sum is odd.
# Therefore it cannot be zero!

# We can see that the Gröbner basis of `S` was computed

@show S.I.gröbner_basis
S.I.algo

# The Gröbner basis is simple to compute in this case as the vector
# of `xi^2 - 1` is already a Gröbner basis.
# However, we still need to divide polynomials by the Gröbner basis
# which can be simplified in this case.

const MP = MultivariatePolynomials
const SS = SemialgebraicSets
struct HypercubeIdeal{V} <: SS.AbstractPolynomialIdeal
    variables::Vector{V}
end
struct HypercubeSet{V} <: SS.AbstractAlgebraicSet
    ideal::HypercubeIdeal{V}
end
MP.variables(set::HypercubeSet) = MP.variables(set.ideal)
MP.variables(ideal::HypercubeIdeal) = ideal.variables
Base.similar(set::HypercubeSet, ::Type) = set
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

# Let's now try to find the correct lower bound:

function min_algebraic_rational(S, d)
    model = SOSModel(solver)
    @variable(model, q, SOSPoly(MP.monomials(x, 0:d)))
    deno = q + 1
    @constraint(model, c, deno * p >= deno, domain = S)
    optimize!(model)
    @show termination_status(model)
end

# With `d = 0`, it's the same as previously

min_algebraic_rational(H, 0)

# But with `d = 1`, we can find the correct lower bound

min_algebraic_rational(H, 1)
