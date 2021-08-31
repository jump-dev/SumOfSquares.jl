# # Sum-of-Squares matrices

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/Sum-of-Squares Matrices.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/Sum-of-Squares Matrices.ipynb)
# **Adapted from**: Examples 3.77 of [BPT12]
#
# [BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
# *Semidefinite Optimization and Convex Algebraic Geometry*.
# Society for Industrial and Applied Mathematics, **2012**.


# ### Introduction

# Consider the symmetric polynomial matrix
# $$P(x) =
# \begin{bmatrix}
#     x^2 - 2x + 2 & x\\
#     x            & x^2
# \end{bmatrix}.$$
# We could like to know whether $P(x)$ is positive semidefinite for all $x \in \mathbb{R}$.

using Test #src
using DynamicPolynomials
@polyvar x
P = [x^2 - 2x + 2 x
            x     x^2]

# A sufficient condition for a symmetric polynomial matrix $P(x)$ to be positive semidefinite for all $x$ is to the existence of a matrix $M(x)$ such that $P(x) = M^\top(x) M(x)$. If such matrix $M$ exists, we say that the matrix is an \emph{sos matrix} (see [Definition 3.76, BPT13]).
# While determining whether $P(x)$ is positive semidefinite for all $x$, is NP-hard (checking nonnegativity of a polynomial is reduced to this problem for $1 \times 1$ matrices), checking whether $P(x)$ is an sos matrix is an sos program.

using SumOfSquares

# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v0.21.6/installation/#Supported-solvers) for a list of the available choices.

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = SOSModel(solver)
mat_cref = @constraint(model, P in PSDCone())
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL #src
termination_status(model) #!src

# While the reformulation of sos matrix to sos polynomial is rather simple, as explained in the "Sum-of-Squares reformulation" section below, there is a technical subtelty about the Newton polytope that if not handled correctly may result in an SDP of large size with bad numerical behavior. For this reason, it is recommended to model sos *matrix* constraints as such as will be shown in this notebook and not do the formulation manually unless there is a specific reason to do so.
#
# As we can verify as follows, only 3 monomials are used using the sos *matrix* constraint.

@test length(certificate_monomials(mat_cref)) == 3 #src
certificate_monomials(mat_cref) #!jl

# ### Sum-of-Squares reformulation
#
# One way to obtain the reduction to an sos program is to create an intermediate vector of variable $y$ and check whether the polynomial $p(x, y) = y^\top P(x) y$ is sos.
# However, special care is required when approximating the Newton polytope of $p(x, y)$.
# Indeed, for instance if the entries of $P(x)$ are quadratic forms then the Newton polytope of $p(x, y)$ is the cartesian product between the Newton polytope of $y^\top y$ and the Newton polytope of $x^\top x$.
# In other words, $p(x, y)$ belongs to a family of quartic forms called biquadratic forms.
# This fact is important when generating the semidefinite program so that only bilinear monomials are used.
# So if the cheap outer approximation is used (instead of the exact polyhedral computation) for the newton polytope then it is important to use a multipartite computation approximation of the newton polytope.
# The multipartie exact approach may perform worse compared to the unipartite exact in certain cases though.
# Consider for instance the polynomial matrix $\mathrm{Diag}(x_1^1, x_2^2)$ for which $p(x, y) = x_1^2y_1^2 + x_2^2y_2^2$.
# For this polynomial, only the monomials $x_1y_1$ and $x_2y_2$ are needed in the SDP reformulation while the multipartite approach,
# as it will compute the Newton polytope as a cartesian product, will not see the dependence between $x$ and $y$ in the presence of monomials and will also select the monomials $x_1y_2$ and $x_2y_1$.

@polyvar y[1:2]
p = vec(y)' * P * vec(y)

# We can see above that `p` is biquadratic polynomial in the variables `x` and `y`.
# Computing the Newton polytope with the cheap outer approximation
# without exploiting this multipartite structure gives the following 6 monomials.

X = monomials(p)
@test Certificate.monomials_half_newton_polytope(X, tuple(), apply_post_filter = false) == [x * y[1], x * y[2], y[1] * y[2], x, y[1], y[2]] #src
Certificate.monomials_half_newton_polytope(X, tuple(), apply_post_filter = false) #!jl

# Exploiting the multipartite structure gives 4 monomials.

@test Certificate.monomials_half_newton_polytope(X, ([x], y), apply_post_filter = false) == [x * y[1], x * y[2], y[1], y[2]] #src
Certificate.monomials_half_newton_polytope(X, ([x], y), apply_post_filter = false) #!jl

# In the example above, there were only 3 monomials, where does the difference come from ?
# Using the monomial basis, the only product of two monomials that is equal to
# `y[2]^2` is `y[2] * y[2]`. As `y[2]^2` is not a monomial of `p`, we can conclude
# that the diagonal entry with row and column corresponding to `y[2]` will be zero
# hence the whole column and row will be zero as well.
# Therefore, we can remove this monomial.

@test Certificate.monomials_half_newton_polytope(X, ([x], y)) == [x * y[1], x * y[2], y[1]] #src
Certificate.monomials_half_newton_polytope(X, ([x], y)) #!jl

# The same reasoning can be used for monomials `y[1]y[2]` and `x` therefore whether
# we exploit the multipartite structure or not, we get only 3 monomials thanks
# to this post filter.

@test Certificate.monomials_half_newton_polytope(X, tuple()) == [x * y[1], x * y[2], y[1]] #src
Certificate.monomials_half_newton_polytope(X, tuple()) #!jl
