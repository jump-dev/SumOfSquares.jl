# # CHSH

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Noncommutative and Hermitian/chsh.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Noncommutative and Hermitian/chsh.ipynb)
# **Contributed by**: Marek Kaluba and Benoît Legat
# **Adapted from**: [Talk](https://jump.dev/assets/jump-dev-workshops/2024/legat.html) at [JuMP-dev 2024](https://jump.dev/meetings/jumpdev2024/)
#
# The goal of this tutorial is to show how to use SumOfSquares.jl with a custom algebra that is **not** defined
# with [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl/) or
# [TypedPolynomials](https://github.com/JuliaAlgebra/TypedPolynomials.jl/).
# Even though some part of SumOfSquares.jl only works for monomials defined by these packages,
# there is an effort to abstract away as much as possible on top of
# [StarAlgebras](https://github.com/JuliaAlgebra/StarAlgebras.jl/).
#
# This illustrate this, in this tutorial, we are doing to define a custom monoid structure implementing the
# rewriting rules of the [CHSH inequality](https://en.wikipedia.org/wiki/CHSH_inequality) using
# [KnuthBendix](https://github.com/kalmarek/KnuthBendix.jl).

using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "monoids.jl"))

import StarAlgebras as SA
import KnuthBendix as KB
import GroupsCore

# The rewriting rule are as follows:

monoid, A, C = trace_monoid(2, 2, A=:A, C=:C)
monoid.rws

# We now define a `StarAlgebra` from [StarAlgebras](https://github.com/JuliaAlgebra/StarAlgebras.jl/)

RM = let monoid = monoid, A = A, C = C, level = 4
    A_l, sizesA = Monoids.wlmetric_ball(A, radius=level)
    C_l, sizesC = Monoids.wlmetric_ball(C, radius=level)
    @time words, sizes = Monoids.wlmetric_ball(
        unique!([a * c for a in A_l for c in C_l]);
        radius=2,
    )
    @info "Sizes of generated balls:" (A, C, combined) =
        (sizesA, sizesC, sizes)
    basis = SA.FixedBasis(words)
    dirac = SA.DiracMStructure(basis, *)
    table = SA.MTable(dirac, (sizes[1], sizes[1]))
    SA.StarAlgebra(monoid, table)
end

# We can convert a monoid element:

A[1], typeof(A[1])

# to an element of the algebra `RM` (so essentially `1 ⋅ A`) as follows:

RM(A[1])

# Then, we can do arithmetic on these algebra element to form the CHSH inequality:

chsh = let A = RM.(A), C = RM.(C)
    A[1] * C[1] + A[1] * C[2] + A[2] * C[1] - A[2] * C[2]
end

# SumOfSquares needs the ambient implicit basis containing all the monoid elements
# so we define it as follows:

struct Full{B} <: SA.ImplicitBasis{B,B} end
Base.in(::B, ::Full{B}) where {B} = true
Base.getindex(::Full{B}, b::B) where {B} = b
import MultivariateBases as MB
MB.implicit_basis(::SA.FixedBasis{B}) where {B} = Full{B}()
function MB.algebra(b::Full{B}) where {B}
    return SA.StarAlgebra(monoid, SA.DiracMStructure(b, *))
end
SA.comparable(::Full) = isless

# The `chsh` polynomial can be rewritten in this basis as follows:

f = SA.AlgebraElement(
    SA.SparseCoefficients(
        [SA.basis(chsh)[k] for (k, _) in SA.nonzero_pairs(SA.coeffs(chsh))],
        [v for (_, v) in SA.nonzero_pairs(SA.coeffs(chsh))],
    ),
    MB.algebra(Full{eltype(SA.basis(chsh))}()),
)

# We pick the SCS solver:

import SCS
scs = optimizer_with_attributes(
    SCS.Optimizer,
    "eps_abs" => 1e-9,
    "eps_rel" => 1e-9,
    "max_iters" => 1000,
)

# We currently need to manually add all SumOfSquares bridges as follows:

model = Model(scs)
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
@variable(model, λ)
@objective(model, Min, λ)

# Adding the SumOfSquares constraint is currently not as user-friendly than
# with monomials:

n = size(SA.mstructure(RM).table, 1)
gram_basis = SA.FixedBasis(SA.basis(chsh).elts[1:n])
cone = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
    SA.basis(chsh),
    [gram_basis],
    [one(f)],
)
@constraint(model, SA.coeffs(λ * one(f) - f, SA.basis(chsh)) in cone);

optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2

# Let's look at the size of the generated SDP:

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
    free -= sum(psd, init=0) do d
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

# We can see that the `FreeBridge` was used because SCS supports free variables and
# affine-in-PSD constraints.

print_active_bridges(model)

# As a workaround for [this MOI issue](https://github.com/jump-dev/MathOptInterface.jl/pull/3001),
# we need to remove the Image bridge

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

# The model is much smaller this time, with 289 equality constraints and 1 free variable.

summary(backend(model).optimizer.model)

# After dualization that is 289 free variables and 1 equality constraints.
# This is a big improvement compared to the 3051 free variables, 18 equality constraints without dualization.

summary(backend(model).optimizer.model.optimizer.dual_problem.dual_model)

# We can see that the bridge used is the `KernelBridge` this time.

print_active_bridges(model)

# We didn't test using the `ImageBridge` with dualization or using the `KernelBridge` without dualization.
# However, these will always produce larger models, which is why the choice o bridge can be need
# automatically once you choose whether you want to dualize the problem or not.
