# # Certificate

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extension/certificate.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extension/certificate.ipynb)
# **Contributed by**: Benoît Legat

# ## Introduction

# Consider the polynomial optimization problem (a variation from [Lasserre2009; Example 2.2](@cite)) of
# minimizing the polynomial $x^3 - x^2 + 2xy - y^2 + y^3$
# over the polyhedron defined by the inequalities $x \ge 0, y \ge 0$ and $x + y \geq 1$.

using Test #src
using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y - y^2 + y^3 + x^3 * y
using SumOfSquares
S = @set x >= 0 && y >= 0 && x^2 + y^2 >= 2

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.
# Note that SumOfSquares generates a *standard form* SDP (i.e., SDP variables
# and equality constraints) while SCS expects a *geometric form* SDP (i.e.,
# free variables and symmetric matrices depending affinely on these variables
# constrained to belong to the PSD cone).
# JuMP will transform the standard from to the geometric form will create the PSD
# variables as free variables and then constrain then to be PSD.
# While this will work, since the dual of a standard from is in in geometric form,
# dualizing the problem will generate a smaller SDP.
# We use therefore `Dualization.dual_optimizer` so that SCS solves the dual problem.

import SCS
using Dualization
solver = dual_optimizer(SCS.Optimizer)

# A Sum-of-Squares certificate that $p \ge \alpha$ over the domain `S`, ensures that $\alpha$ is a lower bound to the polynomial optimization problem.
# The following program searches for the largest lower bound.

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c, p >= α, domain = S)
optimize!(model)

# We can see that the problem is infeasible, meaning that no lower bound was found.

@test termination_status(model) == MOI.INFEASIBLE #src
@test dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE #src
solution_summary(model)

# We now define the Schmüdgen's certificate:

import MultivariateBases as MB
const SOS = SumOfSquares
const SOSC = SOS.Certificate
struct Schmüdgen{IC <: SOSC.AbstractIdealCertificate, CT <: SOS.SOSLikeCone, BT <: SOS.AbstractPolynomialBasis} <: SOSC.AbstractPreorderCertificate
    ideal_certificate::IC
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end

SOSC.cone(certificate::Schmüdgen) = certificate.cone

function SOSC.preprocessed_domain(::Schmüdgen, domain::BasicSemialgebraicSet, p)
    return SOSC.with_variables(domain, p)
end

function SOSC.preorder_indices(::Schmüdgen, domain::SOSC.WithVariables)
    n = length(domain.inner.p)
    if n >= Sys.WORD_SIZE
        error("There are $(2^n - 1) products in Schmüdgen's certificate, they cannot even be indexed with `$Int`.")
    end
    return map(SOSC.PreorderIndex, 1:(2^n-1))
end

function SOSC.multiplier_basis(certificate::Schmüdgen, index::SOSC.PreorderIndex, domain::SOSC.WithVariables)
    q = SOSC.generator(certificate, index, domain)
    return SOSC.maxdegree_gram_basis(certificate.basis, variables(domain), SOSC.multiplier_maxdegree(certificate.maxdegree, q))
end
function SOSC.multiplier_basis_type(::Type{Schmüdgen{IC, CT, BT}}, ::Type{M}) where {IC,CT,BT,M}
    return MB.similar_type(BT, M)
end

function SOSC.generator(::Schmüdgen, index::SOSC.PreorderIndex, domain::SOSC.WithVariables)
    I = [i for i in eachindex(domain.inner.p) if !iszero(index.value & (1 << (i - 1)))]
    return prod([domain.inner.p[i] for i in eachindex(domain.inner.p) if !iszero(index.value & (1 << (i - 1)))])
end

SOSC.ideal_certificate(certificate::Schmüdgen) = certificate.ideal_certificate
SOSC.ideal_certificate(::Type{<:Schmüdgen{IC}}) where {IC} = IC

SOS.matrix_cone_type(::Type{<:Schmüdgen{IC, CT}}) where {IC, CT} = SOS.matrix_cone_type(CT)

# Let's try again with our the Schmüdgen certificate:

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
ideal_certificate = SOSC.Newton(SOSCone(), MB.MonomialBasis, tuple())
certificate = Schmüdgen(ideal_certificate, SOSCone(), MB.MonomialBasis, maxdegree(p))
@constraint(model, c, p >= α, domain = S, certificate = certificate)
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL #src
@test primal_status(model) == MOI.FEASIBLE_POINT #src
@test objective_value(model) ≈ 0.8284 rtol=1e-3 #src
@test length(lagrangian_multipliers(c)) == 7 #src
solution_summary(model)
