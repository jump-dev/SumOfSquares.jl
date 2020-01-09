export sos, chordal_sos, putinar, chordal_putinar

export csp_graph, chordal_csp_graph

function csp_graph(poly::APL, ::FullSpace)
    G = CEG.LabelledGraph{MP.variable_union_type(poly)}()
    for mono in MP.monomials(poly)
        CEG.add_clique!(G, MP.effective_variables(mono))
    end
    return G
end

function csp_graph(poly::APL, domain::AbstractAlgebraicSet)
    G = csp_graph(poly, FullSpace())
    for p in equalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function csp_graph(poly::APL, domain::BasicSemialgebraicSet)
    G = csp_graph(poly, domain.V)
    for p in inequalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function chordal_csp_graph(poly::APL, domain::AbstractBasicSemialgebraicSet)
    return CEG.chordal_extension(csp_graph(poly, domain), CEG.GreedyFillIn())
end

struct Chordal{CT <: SumOfSquares.SOSLikeCone, BT <: PolyJuMP.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end
function get(certificate::MaxDegree, ::GramBasis, poly)
    H, cliques = chordal_csp_graph(poly)
    return map(cliques) do clique
        @assert issorted(clique, rev = true)
        @assert length(Set(clique)) == length(clique)
        return MP.monomials(vars, 0:div(maxdegree, 2))
    end
end

#"""
#    putinar(
#        p::APL,
#        degree::Int,
#        K::AbstractBasicSemialgebraicSet;
#        model::JuMP.Model = SOSModel()
#    )
#
#
#"""
#function putinar(
#        p::APL,
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
#                         p::APL,
#                         degree::Int,
#                         K::AbstractBasicSemialgebraicSet;
#                         model::JuMP.Model = SOSModel()
#                        )
#
#
#"""
#function chordal_putinar(
#                         p::APL,
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
