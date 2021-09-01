# # Sums of Hermitian squares

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/sums_of_hermitian_squares.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/sums_of_hermitian_squares.ipynb)
# **Contributed by**: Benoît Legat

using Test #src

using SumOfSquares

using DynamicPolynomials
@ncpolyvar x y
p = (x + 1.0im * y) * (x - im * y)

import CSDP
model = Model(CSDP.Optimizer)
add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge)
add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge)
add_bridge(model, SumOfSquares.COI.Bridges.Variable.HermitianToSymmetricPSDBridge)
add_bridge(model, SumOfSquares.COI.Bridges.Constraint.SplitZeroBridge)
MOI.Bridges.add_bridge(backend(model).optimizer, PolyJuMP.ZeroPolynomialBridge{Complex{Float64}})
MOI.Bridges.add_bridge(backend(model).optimizer, SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Complex{Float64}})
cone = NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}()
certificate = SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, 2)
c = SumOfSquares.add_constraint(model, p, cone, ideal_certificate = certificate)
optimize!(model)
dec = sos_decomposition(c, 1e-6) #src
@test length(dec.ps) == 1 #src
@test sign(real(first(coefficients(dec.ps[1])))) * dec.ps[1] ≈ x - im * y atol=1e-6 rtol=1e-6 #src
sos_decomposition(c, 1e-6) #!src
