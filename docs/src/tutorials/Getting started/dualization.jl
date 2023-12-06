# # On the importance of Dualization

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/dualization.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/dualization.ipynb)

using Test #src
using DynamicPolynomials
using SumOfSquares

# Sum-of-Squares programs are usually solved by SemiDefinite Programming solvers (SDPs).
# These programs can be represented into two different formats:
# Either the *standard conic form* or *kernel form*:
# ```math
# \begin{aligned}
#   \min\limits_{Q \in \mathbb{S}_n} & \la C, Q \ra\\
#   \text{subject to:} & \la A_i, Q \ra = b_i, \quad i=1,2,\ldots,m\\
#                       & Q \succeq 0,
# \end{aligned}
# ```
# Either the *geometric conic form* or *image form*:
# ```math
# \begin{aligned}
#   \max\limits_{y \in \mathbb{R}^m} & \la b, y \ra\\
#   \text{subject to:} & C \succeq \sum_{i=1}^m A_i y_i\\
#                      & y\ \mathsf{free},
# ```

# In this tutorial, we investigate in which of these two forms a Sum-of-Squares
# constraint should be written into.
# Consider the simple example of trying to see whether the following univariate
# polynomial is Sum-of-Squares:

import SCS
@polyvar x
p = (x + 1)^2 * (x + 2)^2
model_scs = Model(SCS.Optimizer)
con_ref = @constraint(model_scs, p in SOSCone())
optimize!(model_scs)

# As we can see in the log, SCS reports `6` variables and `11` constraints.
# We can also choose to dualize the problem before it is
# passed to SCS as follows:

using Dualization
model_dual_scs = Model(dual_optimizer(SCS.Optimizer))
@objective(model_dual_scs, Max, 0.0)
con_ref = @constraint(model_dual_scs, p in SOSCone())
optimize!(model_dual_scs)

# This time, SCS reports `5` variables and `6` constraints.

print_active_bridges(model_scs)

# This time, we had 

# ## In more details
#
# Consider a polynomial
# ```math
# p(x) = \sum_{\alpha} p_\alpha x^\alpha,
# ```
# a vector of monomials `b(x)` and the set
# ```math
# \mathcal{A}_\alpha = \{\,(\beta, \gamma) \in b(x)^2 \mid x^\beta x^\gamma = x^\alpha\,\}
# ```
# The constraint that there exists a PSD matrix `Q` such that `p(x) = b(x)' * Q * b(x)`
# can be written in standard conic form as follows:
# ```math
# \begin{aligned}
#   \langle \sum_{(\beta, \gamma) \in \mathcal{A}_\alpha} e_\beta e_\gamma^\top, Q \rangle & = p_\alpha, \quad\forall \alpha\\
#   Q & \succeq 0
# \end{aligned}
# ```
# Given an arbitrary choice of elements in each set ``\mathcal{A}_\alpha``:
# ``(\beta_\alpha, \gamma_\alpha) \in \mathcal{A}_\alpha``.
# It can also equivalently be written in the geometric conic form as follows:
# ```math
# \begin{aligned}
#   p_\alpha e_{\beta_\alpha} e_{\gamma_\alpha}^\top +
#   \sum_{(\beta, \gamma) \in \mathcal{A}_\alpha \setminus (\beta_\alpha, \gamma_\alpha)}
#   y_{\beta,\gamma} (e_\beta e_\gamma - e_{\beta_\alpha} e_{\gamma_\alpha}^\top)^\top
#   & \succeq 0\\
#   y_{\beta,\gamma} & \text{ free}
# \end{aligned}
# ```
#
# ## Which one do I choose ?
#
# Let's pick two extreme examples and then try to extrapolate from these.
#
# ### Univariate case 
#
# Suppose `p` is a univariate polynomial of degree $2d$.
# Then `n` will be equal to `d(d + 1)/2` for both the standard and geometric conic forms.
# On the other hand, `m` will be equal to `2d + 1` for the standard conic form and
# `d(d + 1) / 2 - (2d + 1)` for the geometric form case.
# So `m` grows **linearly** for the kernel form but **quadratically** for the image form!
#
# ### Quadratic case 
#
# Suppose `p` is a quadratic form of `d` variables.
# Then `n` will be equal to `d` for both the standard and geometric conic forms.
# On the other hand, `m` will be equal to `d(d + 1)/2` for the standard conic form and
# `0` for the geometric form case.
# So `m` grows **quadratically** for the kernel form but is zero for the image form!
#
# ### In general
#
# In general, if ``s_d`` is the dimension of the space of polynomials of degree `d` then
# ``m = s_{2d}`` for the kernel form and ``m = s_{d}(s_{d} + 1)/2 - s_{2d}`` for the image form.
# As a rule of thumb, the kernel form will have a smaller `m` if `p` has a low number of variables
# and low degree and vice versa.
# Of course, you can always try with and without Dualization and see which one works best.
