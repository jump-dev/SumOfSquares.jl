using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y - y^2 + y^3 + x^3 * y
using SumOfSquares
S = @set x >= 0 && y >= 0 && x^2 + y^2 >= 2

import SCS
using Dualization
solver = dual_optimizer(SCS.Optimizer)

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c, p >= α, domain = S)
optimize!(model)

solution_summary(model)

import StarAlgebras as SA
import MultivariateBases as MB
const SOS = SumOfSquares
const SOSC = SOS.Certificate
struct Schmüdgen{IC<:SOSC.AbstractIdealCertificate,CT<:SOS.SOSLikeCone,B<:SA.AbstractBasis} <: SOSC.AbstractPreorderCertificate
    ideal_certificate::IC
    cone::CT
    basis::B
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
    return MB.explicit_basis_type(BT)
end

function SOSC.generator(::Schmüdgen, index::SOSC.PreorderIndex, domain::SOSC.WithVariables)
    I = [i for i in eachindex(domain.inner.p) if !iszero(index.value & (1 << (i - 1)))]
    return prod([domain.inner.p[i] for i in eachindex(domain.inner.p) if !iszero(index.value & (1 << (i - 1)))])
end

SOSC.ideal_certificate(certificate::Schmüdgen) = certificate.ideal_certificate
SOSC.ideal_certificate(::Type{<:Schmüdgen{IC}}) where {IC} = IC

SOS.matrix_cone_type(::Type{<:Schmüdgen{IC, CT}}) where {IC, CT} = SOS.matrix_cone_type(CT)

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
basis = MB.FullBasis{MB.Monomial,typeof(x * y)}()
ideal_certificate = SOSC.Newton(SOSCone(), basis, tuple())
certificate = Schmüdgen(ideal_certificate, SOSCone(), basis, maxdegree(p))
@constraint(model, c, p >= α, domain = S, certificate = certificate)
optimize!(model)
solution_summary(model)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
