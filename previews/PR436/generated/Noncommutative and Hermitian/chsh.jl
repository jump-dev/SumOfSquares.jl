using SumOfSquares
include(
    joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "monoids.jl"),
)

import StarAlgebras as SA
import KnuthBendix as KB
import GroupsCore

monoid, A, C = trace_monoid(2, 2, A = :A, C = :C)
monoid.rws

RM = let monoid = monoid, A = A, C = C, level = 4
    A_l, sizesA = Monoids.wlmetric_ball(A, radius = level)
    C_l, sizesC = Monoids.wlmetric_ball(C, radius = level)
    @time words, sizes = Monoids.wlmetric_ball(
        unique!([a * c for a in A_l for c in C_l]);
        radius = 2,
    )
    @info "Sizes of generated balls:" (A, C, combined) =
        (sizesA, sizesC, sizes)
    basis = SA.SubBasis(SA.DiracBasis(monoid), words)
    dirac = SA.DiracMStructure(basis, *)
    table = SA.MTable(dirac, (sizes[1], sizes[1]))
    SA.StarAlgebra(monoid, table)
end

A[1], typeof(A[1])

RM(A[1])

chsh = let A = RM.(A), C = RM.(C)
    A[1] * C[1] + A[1] * C[2] + A[2] * C[1] - A[2] * C[2]
end

import MultivariateBases as MB
function MB.algebra(b::SA.DiracBasis{<:Monoids.MonoidElement})
    return SA.StarAlgebra(SA.object(b), b)
end

f = SA.AlgebraElement(
    SA.SparseCoefficients(
        [SA.basis(chsh)[k] for (k, _) in SA.nonzero_pairs(SA.coeffs(chsh))],
        [v for (_, v) in SA.nonzero_pairs(SA.coeffs(chsh))],
    ),
    MB.algebra(parent(SA.basis(chsh))),
)

import SCS
scs = optimizer_with_attributes(
    SCS.Optimizer,
    "eps_abs" => 1e-9,
    "eps_rel" => 1e-9,
    "max_iters" => 1000,
)

model = Model(scs)
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
@variable(model, λ)
@objective(model, Min, λ)

n = size(SA.mstructure(RM).table, 1)
gram_basis = SA.sub_basis(SA.basis(chsh), 1:n)
cone = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
    SA.basis(chsh),
    [gram_basis],
    [one(f)],
)
@constraint(model, SA.coeffs(λ * one(f) - f, SA.basis(chsh)) in cone);

optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2

function _add!(f, psd, model, F, S)
    return append!(
        psd,
        [
            f(MOI.get(model, MOI.ConstraintSet(), ci)) for
            ci in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
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
        return div(d * (d + 1), 2)
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
    println(
        "$free free variables, $(sum(eq, init = 0)) equality constraints, PSD block sizes: $psd",
    )
    return
end
summary(backend(model).optimizer.model)

print_active_bridges(model)

import Dualization
model = Model(Dualization.dual_optimizer(scs))
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
MOI.Bridges.remove_bridge(
    backend(model).optimizer,
    SumOfSquares.Bridges.Constraint.ImageBridge{Float64},
)
@variable(model, λ)
@objective(model, Min, λ)
@constraint(model, SA.coeffs(λ * one(f) - f, SA.basis(chsh)) in cone);
optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2

summary(backend(model).optimizer.model)

summary(backend(model).optimizer.model.optimizer.dual_problem.dual_model)

print_active_bridges(model)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
