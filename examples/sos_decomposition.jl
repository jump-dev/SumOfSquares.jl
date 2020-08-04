using DynamicPolynomials
using SumOfSquares
using JuMP
using SCS

# ------------------------------------------------------------------------------
# A trivial SOS decomposition example
# ------------------------------------------------------------------------------

# The polynomial p = x^2 - x*y^2 + y^4 + 1 is SOS.
# Among infinite others, it has the decomposition:
# p = 3/4*(x-y^2)^2 + 1/4*(x + y)^2 + 1.

# We can use SumOfSquares.jl to find such decompositions.

# First, setup the polynomial of interest. -------------------------------------
@polyvar x y
p = x^2 - x*y^2 + y^4 + 1

# Secondly, constrain the polynomial to be nonnegative. ------------------------
# SumOfSquares.jl transparently reinterprets polyonmial nonnegativity as the
# appropriate SOS certificate for polynomials nonnegative on semialgebraic sets.
model = SOSModel(SCS.Optimizer)
@constraint(model, cref, p >= 0)

# Thirdly, optimize the feasibility problem! -----------------------------------
optimize!(model)

# Lastly, recover a SOS decomposition ------------------------------------------
# In general, SOS decompositions are not unique!
sos_decomposition = SumOfSquares.sos_decomposition(cref, 1e-4)

# Converting, rounding, and simplifying - Huzza, Back where we began!
polynomial(sos_decomposition, Float32)


# ------------------------------------------------------------------------------
# A deeper explanation and the unexplained 1e-4 parameter
# ------------------------------------------------------------------------------

# p = x^2 - x*y^2 + y^4 + 1 can be represented in terms of its Gram matrix as
gram = SumOfSquares.gram_matrix(cref)
gram.basis.monomials' * gram.Q * gram.basis.monomials
# where the matrix gram.Q is positive semidefinite, because p is SOS. If we
# could only get the decomposition gram.Q = V' * V, the SOS decomposition would
# simply be ||V * monomials||^2.

# Unfortunately, we can not use Cholesky decomposition, since gram Q is only
# semidefinite, not definite. Hence, SumOfSquares.jl uses SVD decomposition
# instead and discards small singular values (in our case 1e-4).
