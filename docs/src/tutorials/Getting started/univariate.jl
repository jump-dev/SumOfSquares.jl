# # Mminimization of a univariate polynomial

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/univariate.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/univariate.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
using DynamicPolynomials
using SumOfSquares
import CSDP

# Consider the problem of finding both the minimum value of `p = x^4 - 4x^3 - 2x^2 + 12x + 3` as well as its minimizers.

# We can use SumOfSquares.jl to find such these values as follows.
# We first define the polynomial using DynamicPolynomials.

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 3

# Secondly, we create a Sum-of-Squares program searching for the maximal lower bound `σ` of the polynomial.

model = SOSModel(CSDP.Optimizer)
@variable(model, σ)
@constraint(model, cref, p >= σ)
@objective(model, Max, σ)

# Thirdly, solve the program and find `σ = -6` as lower bound:

optimize!(model)
solution_summary(model)

# We can look at the certificate that `σ = -6` is a lower bound:

sos_dec = sos_decomposition(cref, 1e-4)
expected = x^2 - 2x - 3 #src
@test isapprox(sos_dec.ps, [expected], rtol=1e-4) || isapprox(sos_dec.ps, [-expected], rtol=1e-4) #src

# Indeed, `p + 6 = (x^2 - 2x - 3)^2` so `p ≥ -6`.
#
# ## Extraction of minimizers
#
# We can now find the minimizers from the moment matrix:

ν = moment_matrix(cref)
ν.Q

# This matrix is the convex combination of the moment matrices corresponding to two atomic measures at `-1` and `3`
# which allows us to conclude that `-1` and `3` are global minimizers.

η = extractatoms(ν, 1e-4)
minimizers = [η.atoms[1].center; η.atoms[2].center]

# Below are more details on what we mean by convex combination.
# The moment matrix of the atomic measure at the first minimizer is:

η1 = moment_matrix(dirac(monomials(x, 0:4), x => round(minimizers[1])), ν.basis.monomials)
η1.Q

# The moment matrix of the atomic measure at the second minimizer is:

η2 = moment_matrix(dirac(monomials(x, 0:4), x => round(minimizers[2])), ν.basis.monomials)
η2.Q

# And the moment matrix is the convex combination of both:

Q12 = η1.Q * η.atoms[1].weight + η2.Q * η.atoms[2].weight
@test ν.Q ≈ Q12 rtol=1e-6 #src

# Another way to see this (by linearity of the expectation) is that `ν` is the moment matrix
# of the convex combination of the two atomic measures.

# ## Changing the polynomial basis
#
# The monomial basis used by default can leave a problem quite ill-conditioned for the solver.
# Let's try to use another basis instead:

model = SOSModel(CSDP.Optimizer)
@variable(model, σ)
@constraint(model, cheby_cref, p >= σ, basis = ChebyshevBasisFirstKind)
@objective(model, Max, σ)
optimize!(model)
solution_summary(model)

# Although the gram matrix in the monomial basis:

g = gram_matrix(cref)
@show g.basis
g.Q

# looks different from the gram matrix in the Chebyshev basis:

cheby_g = gram_matrix(cheby_cref)
@show cheby_g.basis
cheby_g.Q

@test polynomial(g) ≈ polynomial(cheby_g) rtol=1e-4 #src

# they both yields the same Sum-of-Squares decomposition:

cheby_sos_dec = sos_decomposition(cheby_cref, 1e-4)
@test isapprox(cheby_sos_dec.ps, [expected], rtol=1e-4) || isapprox(cheby_sos_dec.ps, [-expected], rtol=1e-4) #src
