# # Certificate

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extension/certificate.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extension/certificate.ipynb)
# **Contributed by**: Benoît Legat

# ## Introduction

# Consider the polynomial optimization problem (a variation from [L09, Example 2.2]) of
# minimizing the polynomial $x^3 - x^2 + 2xy - y^2 + y^3$
# over the polyhedron defined by the inequalities $x \ge 0, y \ge 0$ and $x + y \geq 1$.

# [L09] Lasserre, J. B.
# *Moments, positive polynomials and their applications*.
# World Scientific, **2009**.

using Test #src
using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3 + x^3 * y
using SumOfSquares
S = @set x >= 0 && y >= 0 && x^2 + y^2 >= 2

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v0.21.6/installation/#Supported-solvers) for a list of the available choices.

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

# A Sum-of-Squares certificate that $p \ge \alpha$ over the domain `S`, ensures that $\alpha$ is a lower bound to the polynomial optimization problem.
# The following program searches for the largest upper bound and finds zero.

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c, p >= α, domain = S)
optimize!(model)
@show termination_status(model)
@show objective_value(model)

# We now define the Schmüdgen's certificate:

using MultivariateBases
const MB = MultivariateBases
const SOS = SumOfSquares
const SOSC = SOS.Certificate
struct Schmüdgen{IC <: SOSC.AbstractIdealCertificate, CT <: SOS.SOSLikeCone, BT <: SOS.AbstractPolynomialBasis} <: SOSC.AbstractPreorderCertificate
    ideal_certificate::IC
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end

SOSC.get(certificate::Schmüdgen, ::SOSC.Cone) = certificate.cone

function SOSC.get(::Schmüdgen, ::SOSC.PreprocessedDomain, domain::BasicSemialgebraicSet, p)
    return SOSC.DomainWithVariables(domain, variables(p))
end

function SOSC.get(::Schmüdgen, ::SOSC.PreorderIndices, domain::SOSC.DomainWithVariables)
    n = length(domain.domain.p)
    if n >= Sys.WORD_SIZE
        error("There are $(2^n - 1) products in Schmüdgen's certificate, they cannot even be indexed with `$Int`.")
    end
    return map(SOSC.PreorderIndex, 1:(2^n-1))
end

function SOSC.get(certificate::Schmüdgen, ::SOSC.MultiplierBasis, index::SOSC.PreorderIndex, domain::SOSC.DomainWithVariables)
    q = SOSC.get(certificate, SOSC.Generator(), index, domain)
    vars = sort!([domain.variables..., variables(q)...], rev = true)
    unique!(vars)
    return SOSC.maxdegree_gram_basis(certificate.basis, vars, SOSC.multiplier_maxdegree(certificate.maxdegree, q))
end
function SOSC.get(::Type{Schmüdgen{IC, CT, BT}}, ::SOSC.MultiplierBasisType) where {IC, CT, BT}
    return BT
end

function SOSC.get(::Schmüdgen, ::SOSC.Generator, index::SOSC.PreorderIndex, domain::SOSC.DomainWithVariables)
    I = [i for i in eachindex(domain.domain.p) if !iszero(index.value & (1 << (i - 1)))]
    return prod([domain.domain.p[i] for i in eachindex(domain.domain.p) if !iszero(index.value & (1 << (i - 1)))])
end

SOSC.get(certificate::Schmüdgen, ::SOSC.IdealCertificate) = certificate.ideal_certificate
SOSC.get(::Type{<:Schmüdgen{IC}}, ::SOSC.IdealCertificate) where {IC} = IC

SOS.matrix_cone_type(::Type{<:Schmüdgen{IC, CT}}) where {IC, CT} = SOS.matrix_cone_type(CT)

# Let's try again with our the Schmüdgen certificate:

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
ideal_certificate = SOSC.Newton(SOSCone(), MB.MonomialBasis, tuple())
certificate = Schmüdgen(ideal_certificate, SOSCone(), MB.MonomialBasis, maxdegree(p))
@constraint(model, c, p >= α, domain = S, certificate = certificate)
optimize!(model)
@show termination_status(model)
@show objective_value(model)
@test length(lagrangian_multipliers(c)) == 7 #src
