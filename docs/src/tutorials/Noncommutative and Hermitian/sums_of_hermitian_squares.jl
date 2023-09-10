# # Sums of Hermitian squares

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Noncommutative and Hermitian/sums_of_hermitian_squares.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Noncommutative and Hermitian/sums_of_hermitian_squares.ipynb)
# **Contributed by**: Benoît Legat

using Test #src

using SumOfSquares

using DynamicPolynomials
@ncpolyvar x y
p = (x + im * y) * (x - im * y)

# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.8/installation/#Supported-solvers) for a list of the available choices.

import Clarabel
model = Model(Clarabel.Optimizer)
cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
con_ref = @constraint(model, p in cone)
optimize!(model)
dec = sos_decomposition(con_ref, 1e-6) #src
@test length(dec.ps) == 1 #src
@test dec.ps[1]' * dec.ps[1] ≈ p atol=1e-6 rtol=1e-6 #src
@test sign(real(first(coefficients(dec.ps[1])))) * dec.ps[1] ≈ y + im * x atol=1e-6 rtol=1e-6 #src
sos_decomposition(con_ref, 1e-6)
