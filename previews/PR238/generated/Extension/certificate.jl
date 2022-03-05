using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3 + x^3 * y
using SumOfSquares
S = @set x >= 0 && y >= 0 && x^2 + y^2 >= 2

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c, p >= α, domain = S)
optimize!(model)
@show termination_status(model)
@show objective_value(model)

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

SOSC.cone(certificate::Schmüdgen) = certificate.cone

function SOSC.preprocessed_domain(::Schmüdgen, domain::BasicSemialgebraicSet, p)
    return SOSC.DomainWithVariables(domain, variables(p))
end

function SOSC.preorder_indices(::Schmüdgen, domain::SOSC.DomainWithVariables)
    n = length(domain.domain.p)
    if n >= Sys.WORD_SIZE
        error("There are $(2^n - 1) products in Schmüdgen's certificate, they cannot even be indexed with `$Int`.")
    end
    return map(SOSC.PreorderIndex, 1:(2^n-1))
end

function SOSC.multiplier_basis(certificate::Schmüdgen, index::SOSC.PreorderIndex, domain::SOSC.DomainWithVariables)
    q = SOSC.generator(certificate, index, domain)
    vars = sort!([domain.variables..., variables(q)...], rev = true)
    unique!(vars)
    return SOSC.maxdegree_gram_basis(certificate.basis, vars, SOSC.multiplier_maxdegree(certificate.maxdegree, q))
end
function SOSC.multiplier_basis_type(::Type{Schmüdgen{IC, CT, BT}}) where {IC, CT, BT}
    return BT
end

function SOSC.generator(::Schmüdgen, index::SOSC.PreorderIndex, domain::SOSC.DomainWithVariables)
    I = [i for i in eachindex(domain.domain.p) if !iszero(index.value & (1 << (i - 1)))]
    return prod([domain.domain.p[i] for i in eachindex(domain.domain.p) if !iszero(index.value & (1 << (i - 1)))])
end

SOSC.ideal_certificate(certificate::Schmüdgen) = certificate.ideal_certificate
SOSC.ideal_certificate(::Type{<:Schmüdgen{IC}}) where {IC} = IC

SOS.matrix_cone_type(::Type{<:Schmüdgen{IC, CT}}) where {IC, CT} = SOS.matrix_cone_type(CT)

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
ideal_certificate = SOSC.Newton(SOSCone(), MB.MonomialBasis, tuple())
certificate = Schmüdgen(ideal_certificate, SOSCone(), MB.MonomialBasis, maxdegree(p))
@constraint(model, c, p >= α, domain = S, certificate = certificate)
optimize!(model)
@show termination_status(model)
@show objective_value(model)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

