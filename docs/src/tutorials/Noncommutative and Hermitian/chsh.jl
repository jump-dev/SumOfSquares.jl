import QuantumStuff: trace_monoid, Monoids
using StarAlgebras

M, A, C = trace_monoid(2, 2, A=:A, C=:C)

RM = let M = M, A = A, C = C, level = 4
    A_l, sizesA = Monoids.wlmetric_ball(A, radius=level)
    C_l, sizesC = Monoids.wlmetric_ball(C, radius=level)

    # starAlg(M, 1, half = unique!([a*c for a in A_l for c in C_l]))

    @time words, sizes = Monoids.wlmetric_ball(
        unique!([a * c for a in A_l for c in C_l]);
        radius=2,
    )
    @info "Sizes of generated balls:" (A, C, combined) = (sizesA, sizesC, sizes)

    b = @time StarAlgebras.FixedBasis(words, StarAlgebras.DiracMStructure(*), (UInt32(sizes[1]), UInt32(sizes[1])))
    StarAlgebra(M, b)
end

A = RM.(A)
C = RM.(C)
chsh = A[1] * C[1] + A[1] * C[2] + A[2] * C[1] - A[2] * C[2]

import StarAlgebras as SA
struct Full{B} <: SA.ImplicitBasis{B,B} end
Base.getindex(::Full{B}, b::B) where {B} = b
import MultivariateBases as MB
MB.implicit_basis(::SA.FixedBasis{B}) where {B} = Full{B}()
MB.algebra(b::Full{B}) where {B} = SA.StarAlgebra(M, b)
SA.mstructure(::Full) = SA.DiracMStructure(*)

b = basis(chsh)
import StarAlgebras as SA
f = SA.AlgebraElement(
    SA.SparseCoefficients(
        [b[k] for (k, _) in SA.nonzero_pairs(coeffs(chsh))],
        [v for (_, v) in SA.nonzero_pairs(coeffs(chsh))],
    ),
    SA.StarAlgebra(
        parent(chsh).object,
        Full{eltype(b)}()
    ),
)
n = size(b.table, 1)
gram_basis = @time StarAlgebras.FixedBasis(b.elts[1:n], StarAlgebras.DiracMStructure(*));
one(f)
SA.coeffs(f, b)
using SumOfSquares
function SumOfSquares._term_element(a, mono::Monoids.MonoidElement)
    SA.AlgebraElement(
        SA.SparseCoefficients((mono,), (a,)),
        MB.algebra(Full{typeof(mono)}()),
    )
end

cone = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
    b,
    [gram_basis],
    [one(f)],
)
import SCS
scs = optimizer_with_attributes(
    SCS.Optimizer,
    "eps_abs" => 1e-9,
    "eps_rel" => 1e-9,
    "max_iters" => 1000,
)

import Dualization
#model = Model(Dualization.dual_optimizer(scs))
model = Model(scs)
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
MOI.Bridges.remove_bridge(backend(model).optimizer, SumOfSquares.Bridges.Constraint.ImageBridge{Float64})
@variable(model, λ)
@objective(model, Min, λ)
@constraint(model, SA.coeffs(λ * one(f) - f, b) in cone);
optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2
function _add!(f, psd, model, F, S)
    append!(psd, [
        f(MOI.get(model, MOI.ConstraintSet(), ci))
        for ci in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
    ],
    )
end
function summary(model)
    free = MOI.get(model, MOI.NumberOfVariables())
    psd = Int[]
    F = MOI.VectorOfVariables
    S = MOI.PositiveSemidefiniteConeTriangle
    _add!(MOI.side_dimension, psd, model, F, S)
    S = SCS.ScaledPSDCone
    _add!(MOI.side_dimension, psd, model, F, S)
    free -= sum(psd, init = 0) do d
        div(d * (d + 1), 2)
    end
    F = MOI.VectorAffineFunction{Float64}
    S = MOI.PositiveSemidefiniteConeTriangle
    _add!(MOI.side_dimension, psd, model, F, S)
    S = SCS.ScaledPSDCone
    _add!(MOI.side_dimension, psd, model, F, S)
    eq = Int[]
    F = MOI.VectorAffineFunction{Float64}
    S = MOI.Zeros
    _add!(MOI.dimension, eq, model, F, S)
    F = MOI.ScalarAffineFunction{Float64}
    S = MOI.EqualTo{Float64}
    _add!(MOI.dimension, eq, model, F, S)
    return free, psd, sum(eq, init = 0)
end
summary(backend(model).optimizer.model)
summary(backend(model).optimizer.model.optimizer.dual_problem.dual_model.model)
print_active_bridges(model)
