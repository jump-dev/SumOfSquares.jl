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
    return CEG.chordal_extension(csp_graph(poly, domain), CEG.GreedyFillIn())
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
        @assert issorted(clique, rev = true)
        @assert length(Set(clique)) == length(clique)
        return MP.monomials(clique, 0:div(certificate.maxdegree, 2))
    end
end

#"""
#    putinar(
#        p::MP.APL,
#        degree::Int,
#        K::AbstractBasicSemialgebraicSet;
#        model::JuMP.Model = SOSModel()
#    )
#
#
#"""
#function putinar(
#        p::MP.APL,
#        degree::Int,
#        K::AbstractBasicSemialgebraicSet;
#        model::JuMP.Model = SOSModel()
#    )
#
#    if K isa FullSpace
#        println("Unbounded domain. Ignoring degree = $degree .")
#        model = sos(p, model)
#    else
#        vars = variables(p)
#        if K isa BasicSemialgebraicSet
#            for ineq in inequalities(K)
#                append!(vars, variables(ineq))
#            end
#        end
#        for eq in equalities(K)
#            append!(vars, variables(eq))
#        end
#        vars = Tuple(unique!(sort!(vars, rev = true)))
#
#        while !isempty(equalities(K))
#            eq = pop!(equalities(K))
#            deg = degree - MP.maxdegree(eq)
#            if deg >= 0
#                 mvec = MP.monomials(vars, 0:deg)
#                 ppi = @variable model variable_type=Poly(mvec)
#                 println(ppi)
#                 p = p - ppi*eq
#            end
#        end
#        @constraint model p in SOSCone() domain = K maxdegree = degree
#    end
#    return model
#end
#
#
#"""
#    chordal_putinar(
#                         p::MP.APL,
#                         degree::Int,
#                         K::AbstractBasicSemialgebraicSet;
#                         model::JuMP.Model = SOSModel()
#                        )
#
#
#"""
#function chordal_putinar(
#                         p::MP.MP.APL,
#                         degree::Int,
#                         K::AbstractBasicSemialgebraicSet;
#                         model::JuMP.Model = SOSModel()
#                        )
#
#    if K isa FullSpace
#        println("Unbounded domain. Ignoring degree = $degree .")
#        model = chordal_sos(p, model)
#    else
#        H, cliques = chordal_csp_graph(p, K)
#
#        for clique in cliques
#
#            vars = Tuple(unique!(sort!(clique, rev = true)))
#            mvec = MP.monomials(vars, 0:degree)
#            pp = @variable model variable_type=Poly(mvec)
#
#            p = p - pp
#
#            dummy = sum(rand()*mvec[i] for i=1:length(mvec))
#            Ki = BasicSemialgebraicSet{Float64, typeof(dummy)}()
#
#            if K isa BasicSemialgebraicSet
#                for ineq in inequalities(K)
#                    if variables(ineq) ⊆ clique
#                        addinequality!(Ki, ineq)
#                    end
#                end
#            end
#
#            for eq in equalities(K)
#                if variables(eq) ⊆ clique
#                    deg = degree - MP.maxdegree(eq)
#                    if deg >= 0
#                        mvec = MP.monomials(vars, 0:deg)
#                        ppi = @variable model variable_type=Poly(mvec)
#                        pp = pp - ppi*eq
#                    end
#                end
#            end
#            @constraint model pp in SOSCone() domain = Ki maxdegree = degree
#        end
#        @constraint model p == 0
#    end
#    return model
#end
#
