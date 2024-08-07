# # Univariate Solver

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extension/univariate_solver.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extension/univariate_solver.ipynb)
# **Contributed by**: Benoît Legat

# When using an SDP solver, the Sum-of-Squares constraint is bridged
# into a semidefinite constraint. This reformulation is done only if the
# solver does not support the Sum-of-Squares constraint. We show in this tutorial how to define a solver that supports such constraint.
# The *same* model with be then solved with or without reformulation depending on the solver.

# We consider a solver that would support finding the SOS decomposition
# of nonnegative univariate polynomials.

using Test #src

module MyUnivariateSolver

import LinearAlgebra
import MathOptInterface as MOI
import MultivariatePolynomials as MP
import MultivariateBases as MB
import SumOfSquares as SOS

function decompose(p::MP.AbstractPolynomial, tol=1e-6)
    vars = MP.effective_variables(p)
    if length(vars) != 1
        error("$p is not univariate")
    end
    x = first(vars)
    lead = MP.leading_coefficient(p)
    if !isone(lead)
        p = p / lead
    end
    deg = MP.maxdegree(p)
    if isodd(deg)
        return
    end
    d = div(deg, 2)
    companion = zeros(2d, 2d)
    for i in 0:(2d-1)
        if i > 0
            companion[i + 1, i] = 1.0
        end
        companion[i + 1, end] = -MP.coefficient(p, x^i)
    end
    F = LinearAlgebra.schur(complex(companion))
    q = one(p)
    i = 1
    while i <= length(F.values)
        root = F.values[i]
        q *= (x - root)
        if !isapprox(real(root), real(F.values[i+1]), rtol=tol, atol=tol)
            return # Cannot happen for complex conjugate root so it means that we have a root which does not have an even multiplicity This means that the polynomial is not nonnegative
        end
        i += 2
    end
    q1 = MP.map_coefficients(real, q)
    q2 = MP.map_coefficients(imag, q)
    return SOS.SOSDecomposition([MB.algebra_element(q1), MB.algebra_element(q2)])
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    p::Union{Nothing,MP.AbstractPolynomial}
    decomposition::Union{Nothing,SOS.SOSDecomposition}
    tol::Float64
    function Optimizer()
        return new(nothing, nothing, 1e-6)
    end
end

MOI.is_empty(optimizer::Optimizer) = optimizer.p === nothing
function MOI.empty!(optimizer::Optimizer)
    optimizer.p = nothing
    return
end

function MOI.supports_constraint(::Optimizer, ::Type{<:MOI.VectorAffineFunction}, ::Type{<:SOS.SOSPolynomialSet{SOS.FullSpace}})
    return true
end
function MOI.add_constraint(optimizer::Optimizer, func::MOI.VectorAffineFunction, set::SOS.SOSPolynomialSet{SOS.FullSpace})
    if optimizer.p !== nothing
        error("Only one constraint is supported")
    end
    if !isempty(func.terms)
        error("Only supports constant polynomials")
    end
    optimizer.p = MP.polynomial(MB.algebra_element(func.constants, set.basis))
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(1) # There will be only ever one constraint so the index does not matter.
end

MOI.supports_incremental_interface(::Optimizer) = true
function MOI.copy_to(optimizer::Optimizer, model::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(optimizer, model)
end
function MOI.optimize!(optimizer::Optimizer)
    optimizer.decomposition = decompose(optimizer.p, optimizer.tol)
end

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    if optimizer.decomposition === nothing
        return MOI.INFEASIBLE
    else
        return MOI.OPTIMAL
    end
end

function MOI.get(optimizer::Optimizer, ::MOI.PrimalStatus)
    if optimizer.decomposition === nothing
        return MOI.NO_SOLUTION
    else
        return MOI.FEASIBLE_POINT
    end
end

function MOI.get(optimizer::Optimizer, ::SOS.SOSDecompositionAttribute, ::MOI.ConstraintIndex)
    return optimizer.decomposition
end

end

# We define the same function both for our new solver and SDP solvers.

using SumOfSquares
function decompose(p, solver)
    model = Model(solver)
    con = @constraint(model, p in SOSCone())
    optimize!(model)
    @assert primal_status(model) == MOI.FEASIBLE_POINT
    return sos_decomposition(con)
end

# We consider the following univariate polynomial from
# [Blekherman2012; Example 3.35](@cite).

# Using our solver we find the following decomposition in two squares.

using DynamicPolynomials
@polyvar x
p = x^4 + 4x^3 + 6x^2 + 4x + 5
dec = decompose(p, MyUnivariateSolver.Optimizer)

# We can verify as follows that it gives a correct decomposition:

@test polynomial(dec) ≈ p #src
polynomial(dec)

# We can also use a semidefinite solver:

import CSDP
dec = decompose(p, CSDP.Optimizer)
@test length(dec.ps) == 3 #src
@test polynomial(dec) ≈ p #src

# The decomposition is different, it is the sum of 3 squares.
# However, it is also valid:

polynomial(dec)
