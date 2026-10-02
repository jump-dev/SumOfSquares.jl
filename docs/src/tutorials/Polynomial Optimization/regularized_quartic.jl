# # Quartically regularized polynomials

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/regularized_quartic.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/regularized_quartic.ipynb)
# **Adapted from**: Section 3.3.1 of [Cartis2026](@cite)

# Higher-order optimization methods minimize, at each iteration, a cubic
# Taylor model regularized by a quartic term:
# ```math
# m(s) = f_0 + g^\top s + \frac{1}{2} H[s]^2 + \frac{1}{6} T[s]^3 + \frac{\sigma}{4} \|s\|_2^4.
# ```
# For any $\sigma > 0$, this model is bounded below so $m(s) - m^*$ is
# nonnegative where $m^*$ is its global minimum. If $m(s) - m^*$ is moreover
# a sum of squares then $m^*$ can be computed exactly with a single
# Sum-of-Squares program, i.e., the zeroth level of the Lasserre hierarchy.
# [Cartis2026](@cite) shows that this is the case for $\sigma$ large enough
# and that, when $T = 0$, it is the case for all $\sigma > 0$.
# They also show that the geometry of the regularization matters:
# with the *separable* quartic regularization
# $\|s\|_4^4 = s_1^4 + \cdots + s_n^4$ instead of the Euclidean one
# $\|s\|_2^4 = (s_1^2 + \cdots + s_n^2)^2$,
# $m(s) - m^*$ may fail to be a sum of squares for every $\sigma > 0$.
# In this tutorial, we verify this numerically on their example.

using Test #src
using LinearAlgebra
using DynamicPolynomials
@polyvar s[1:3]

# We consider the following nonconvex quadratic part whose Hessian has
# eigenvalues $-4$, $-4$ and $20$:

quad = 2 * sum(s .^ 2) + 8 * (s[1] * s[2] + s[1] * s[3] + s[2] * s[3])

# We need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) for a list of the available choices.

using SumOfSquares
import Clarabel
solver = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

# The following function computes the largest $\gamma$ such that
# `multiplier * (p - γ)` is a sum of squares. With the default `multiplier`,
# this is the zeroth level of the Lasserre hierarchy. It returns
# the lower bound $\gamma$ and the constraint reference, which we will use to
# extract minimizers from the moment matrix.

function sos_lower_bound(p; multiplier = 1)
    model = SOSModel(solver)
    @variable(model, γ)
    @objective(model, Max, γ)
    con_ref = @constraint(model, multiplier * (p - γ) >= 0)
    optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL
    return value(γ), con_ref
end

# ## Separable quartic regularization

# We first consider the separable regularization with $\sigma = 4$.
# The resulting polynomial was shown in [Ahmadi2023; Theorem 3.3](@cite)
# to be such that $m_A(s) - m_A^*$ is not a sum of squares.

m_A = quad + sum(s .^ 4)

# The zeroth level of the Lasserre hierarchy gives the following lower bound:

γ_A, _ = sos_lower_bound(m_A)
@test γ_A ≈ -2.38595 rtol = 1e-4 #src
γ_A

# Is this lower bound tight ? Recall that a polynomial $p$ is nonnegative if
# $(s_1^2 + s_2^2 + s_3^2) p(s)$ is a sum of squares, see the
# [Motzkin](@ref) tutorial. Using this multiplier, we get a larger
# lower bound:

γ_A_mult, con_ref = sos_lower_bound(m_A, multiplier = sum(s .^ 2))
@test γ_A_mult ≈ -2.14428 rtol = 1e-4 #src
γ_A_mult

# We can extract the minimizers from the moment matrix:

atoms = atomic_measure(moment_matrix(con_ref), 1e-4)

# We recover the six global minimizers $\pm(-1.1498, 0.6674, 0.6674)$ and
# their permutations. The moment matrix also contains an atom at the origin
# with zero weight. It comes from the multiplier which vanishes at the
# origin, so we ignore it.

minimizers = [atom.center for atom in atoms.atoms if atom.weight > 1e-4]
@test length(minimizers) == 6 #src

# Evaluating $m_A$ at the minimizers confirms that this is the global minimum
# and that the zeroth level of the hierarchy has a gap:

m_A_min = minimum(x -> m_A(s => x), minimizers)
@test m_A_min ≈ γ_A_mult rtol = 1e-4 #src
m_A_min

# As shown in [Cartis2026; Lemma 3.3](@cite), the gap does not close as
# $\sigma$ increases. Indeed, with $\tilde{s} = \sqrt{\sigma}s/2$, we have
# $\mathrm{quad}(s) + \sigma \|s\|_4^4/4 = 4 m_A(\tilde{s})/\sigma$
# so the ratio between the lower bound and the minimum is independent of
# $\sigma$:

for σ in [1, 10, 100]
    m = quad + σ / 4 * sum(s .^ 4)
    γ, _ = sos_lower_bound(m)
    @test γ * σ / 4 ≈ γ_A rtol = 1e-4 #src
    println("σ = $σ: lower bound = $γ, minimum = $(4m_A_min / σ)")
end

# ## Euclidean quartic regularization

# Let us now replace the separable regularization with the Euclidean one:

σ = 4
m_E = quad + σ / 4 * sum(s .^ 2)^2

# This is a quadratic model with Euclidean quartic regularization ($T = 0$) so
# $m_E(s) - m_E^*$ is a sum of squares for every $\sigma > 0$.
# Indeed, it can be checked that
# ```math
# m_E(s) + \frac{4}{\sigma}
# = \frac{\sigma}{4}\left(\|s\|_2^2 - \frac{4}{\sigma}\right)^2 + 4(s_1 + s_2 + s_3)^2.
# ```

@test m_E + 4 / σ ≈ σ / 4 * (sum(s .^ 2) - 4 / σ)^2 + 4 * sum(s)^2 #src

# The zeroth level of the hierarchy finds this lower bound $-4/\sigma = -1$:

γ_E, con_ref = sos_lower_bound(m_E)
@test γ_E ≈ -4 / σ rtol = 1e-6 #src
γ_E

# The global minimizers form the circle $s_1 + s_2 + s_3 = 0$,
# $\|s\|_2 = 2/\sqrt{\sigma}$. Even if there are infinitely many of them,
# the low rank decomposition of the moment matrix still finds a few points
# of this circle:

atoms = atomic_measure(moment_matrix(con_ref), 1e-4)

# We can verify that they belong to the circle:

for atom in atoms.atoms
    x = atom.center
    @test sum(x) ≈ 0 atol = 1e-4 #src
    @test norm(x) ≈ 2 / √σ rtol = 1e-4 #src
    println("sum = $(sum(x)), norm = $(norm(x))")
end

# ## Adding a cubic term

# We now add a small cubic term $0.1(s_1^3 + s_2^3 + s_3^3)$ to both models.
# With the separable regularization, there is still a gap between the
# lower bound of the zeroth level and the global minimum:

cubic = 0.1 * sum(s .^ 3)
m_A2 = m_A + cubic
γ_A2, _ = sos_lower_bound(m_A2)
@test γ_A2 ≈ -2.39899 rtol = 1e-4 #src
γ_A2

# The multiplier closes the gap again and gives the global minimum:

γ_A2_mult, con_ref = sos_lower_bound(m_A2, multiplier = sum(s .^ 2))
@test γ_A2_mult ≈ -2.24095 rtol = 1e-4 #src
γ_A2_mult

# With the Euclidean regularization, the model is now a cubic model so
# [Cartis2026; Theorem 2.3](@cite) only guarantees that $m(s) - m^*$ is a sum of
# squares for $\sigma$ large enough. We compare the lower bound found by
# the zeroth level of the hierarchy with the value of $m$ at the minimizers
# extracted from the moment matrix. If they match, then the lower bound is the
# global minimum.

for σ in [0.1, 1, 4, 10, 40]
    m = quad + cubic + σ / 4 * sum(s .^ 2)^2
    γ, c = sos_lower_bound(m)
    μ = atomic_measure(moment_matrix(c), 1e-4)
    m_min = minimum(atom -> m(s => atom.center), μ.atoms)
    @test γ ≈ m_min rtol = 1e-4 #src
    println("σ = $σ: lower bound = $γ, value at extracted minimizers = $m_min")
end

# For this example, the zeroth level is exact even for small values of
# $\sigma$. The theoretical threshold of [Cartis2026; Theorem 2.3](@cite) is
# sufficient but not necessary.
