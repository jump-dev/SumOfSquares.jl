export sos, chordal_sos, putinar, chordal_putinar

export csp_graph, chordal_csp_graph

const CEG = ChordalExtensionGraph

function csp_graph(poly::MP.APL, ::FullSpace)
    G = CEG.LabelledGraph{MP.variable_union_type(poly)}()
    for mono in MP.monomials(poly)
        CEG.add_clique!(G, MP.effective_variables(mono))
    end
    return G
end

function csp_graph(poly::MP.APL, domain::AbstractAlgebraicSet)
    G = csp_graph(poly, FullSpace())
    for p in equalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function csp_graph(poly::MP.APL, domain::BasicSemialgebraicSet)
    G = csp_graph(poly, domain.V)
    for p in inequalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function chordal_csp_graph(poly::MP.APL, domain::AbstractBasicSemialgebraicSet)
    H, cliques = CEG.chordal_extension(csp_graph(poly, domain), CEG.GreedyFillIn())
    for clique in cliques
        sort!(clique, rev=true)
        unique!(clique)
    end
    return H, cliques

end

struct ChordalPutinar{CT <: SumOfSquares.SOSLikeCone, BT <: PolyJuMP.AbstractPolynomialBasis} <: AbstractPreorderCertificate
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end

struct ChordalDomain{S, V}
    domain::S
    cliques::Vector{V}
end

function get(::ChordalPutinar, ::PreprocessedDomain, domain::BasicSemialgebraicSet, p)
    H, cliques = chordal_csp_graph(p, domain)
    return ChordalDomain(domain, cliques)
end

function get(::ChordalPutinar, index::PreorderIndices, domain::ChordalDomain)
    return map(PreorderIndex, eachindex(domain.domain.p))
end

function get(certificate::ChordalPutinar, ::MultiplierBasis, index::PreorderIndex, domain::ChordalDomain)
    q = domain.domain.p[index.value]
    maxdegree_s2 = certificate.maxdegree - MP.maxdegree(q)
    # If maxdegree_s2 is odd, `div(maxdegree_s2, 2)` would make s^2 have degree up to maxdegree_s2-1
    # for this reason, we take `div(maxdegree_s2 + 1, 2)` so that s^2 have degree up to maxdegree_s2+1
    maxdegree_s = div(maxdegree_s2 + 1, 2)
    return [MP.monomials(clique, 0:maxdegree_s) for clique in domain.cliques if variables(q) ⊆ clique]
end

function get(::ChordalPutinar, ::Generator, index::PreorderIndex, domain::ChordalDomain)
    return domain.domain.p[index.value]
end

get(certificate::ChordalPutinar, ::IdealCertificate) = ChordalIdeal(certificate.cone, certificate.basis, certificate.maxdegree)
get(::Type{<:ChordalPutinar{CT, BT}}, ::IdealCertificate) where {CT, BT} = ChordalIdeal{CT, BT}

SumOfSquares.matrix_cone_type(::Type{<:ChordalPutinar{CT}}) where {CT} = SumOfSquares.matrix_cone_type(CT)

struct ChordalIdeal{CT <: SumOfSquares.SOSLikeCone, BT <: PolyJuMP.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end
function get(certificate::ChordalIdeal, ::GramBasis, poly)
    H, cliques = chordal_csp_graph(poly, FullSpace())
    return map(cliques) do clique
        return MP.monomials(clique, 0:div(certificate.maxdegree, 2))
    end
end
