# # Noncommutative variables

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Sum-of-Squares Matrices.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Sum-of-Squares Matrices.ipynb)
# **Adapted from**: Examples 2.11 and 2.2 of [BKP16]
#
# [BKP16] Sabine Burgdorf, Igor Klep, and Janez Povh.
# *Optimization of polynomials in non-commuting variables*.
# Berlin: Springer, 2016.

# ## Example 2.11
#
# We consider the Example 2.11 of [BKP16] in which the polynomial with noncommutative variables
# $(x * y + x^2)^2 = x^4 + x^3y + xyx^2 + xyxy$ is tested to be sum-of-squares.

using DynamicPolynomials
@ncpolyvar x y
p = (x * y + x^2)^2

import CSDP
using JuMP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
using SumOfSquares
model = Model(optimizer_constructor)
con_ref = @constraint(model, p in SOSCone())

optimize!(model)

# We see that both the monomials `xy` and `yx` are considered separately, this is a difference with the commutative version.

certificate_basis(con_ref)

# We see that the solution correctly uses the monomial `xy` instead of `yx`. We also identify that only the monomials `x^2` and `xy` would be needed. This would be dectected by the Newton chip method of [Section 2.3, BKP16].

gram_matrix(con_ref).Q

# When asking for the SOS decomposition, the numerically small entries makes the solution less readable.

sos_decomposition(con_ref)

# They are however easily discarded by using a nonzero tolerance:

sos_decomposition(con_ref, 1e-6)

# ## Example 2.2
#
# We consider now the Example 2.2 of [BKP16] in which the polynomial with noncommutative variables
# $(x + x^{10}y^{20}x^{10})^2$ is tested to be sum-of-squares.

using DynamicPolynomials
@ncpolyvar x y
n = 10
p = (x + x^n * y^(2n) * x^n)^2

using SumOfSquares
model = Model(optimizer_constructor)
con_ref = @constraint(model, p in SOSCone())

optimize!(model)

# Only two monomials were considered for the basis of the gram matrix thanks to the Augmented Newton chip method detailed in [Section 2.4, BKP16].

certificate_basis(con_ref)

gram_matrix(con_ref).Q

sos_decomposition(con_ref, 1e-6)
