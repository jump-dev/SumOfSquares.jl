using Test #src

using SumOfSquares

import ComplexOptInterface
const COI = ComplexOptInterface
function SumOfSquares.matrix_constructor(::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, T::Type)
    return (Q, basis) -> SumOfSquares.MultivariateMoments.VectorizedHermitianMatrix{eltype(Q), T}(Q, basis)
end
function SumOfSquares.matrix_constructor(::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, ::Type{Complex{T}}) where T
    return SumOfSquares.matrix_constructor(COI.HermitianPositiveSemidefiniteConeTriangle, T)
end
function SumOfSquares.matrix_cone(::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, d)
    return COI.HermitianPositiveSemidefiniteConeTriangle(d)
end

using DynamicPolynomials
@ncpolyvar x y
p = (x + 1.0im * y) * (x - im * y)

import CSDP

model = Model(CSDP.Optimizer)
add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge)
add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge)
add_bridge(model, COI.Bridges.Variable.HermitianToSymmetricPSDBridge)
add_bridge(model, COI.Bridges.Constraint.SplitZeroBridge)
MOI.Bridges.add_bridge(backend(model).optimizer, PolyJuMP.ZeroPolynomialBridge{Complex{Float64}})
MOI.Bridges.add_bridge(backend(model).optimizer, SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Complex{Float64}})
cone = NonnegPolyInnerCone{COI.HermitianPositiveSemidefiniteConeTriangle}()
certificate = SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, 2)
c = SumOfSquares.add_constraint(model, p, cone, ideal_certificate = certificate)
optimize!(model)
dec = sos_decomposition(c, 1e-6) #src
@test length(dec.ps[1]) == 1 #src
@test sign(real(first(coefficients(dec.ps[1])))) * dec.ps[1] ≈ x - im * y atol=1e-6 rtol=1e-6 #src
sos_decomposition(c, 1e-6) #!src
