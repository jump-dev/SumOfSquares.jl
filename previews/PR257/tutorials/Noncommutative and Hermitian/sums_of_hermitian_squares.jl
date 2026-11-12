# # Sums of Hermitian squares

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Noncommutative and Hermitian/sums_of_hermitian_squares.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Noncommutative and Hermitian/sums_of_hermitian_squares.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
import StarAlgebras as SA #src

using SumOfSquares

using DynamicPolynomials
@ncpolyvar x y
p = (x + im * y) * (x - im * y)

import CSDP
model = Model(CSDP.Optimizer)
cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
con_ref = @constraint(model, p in cone)
optimize!(model)
dec = sos_decomposition(con_ref, 1e-6) #src
@test length(dec.ps) == 1 #src
@test polynomial(dec.ps[1])' * polynomial(dec.ps[1]) ≈ p atol=1e-6 rtol=1e-6 #src
@test sign(real(first(SA.coeffs(dec.ps[1])))) * dec.ps[1] ≈ y + im * x atol=1e-6 rtol=1e-6 #src
sos_decomposition(con_ref, 1e-6)
