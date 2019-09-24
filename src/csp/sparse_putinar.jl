export chordal_sos, chordal_putinar

"""
    chordal_sos(p::APL; model::JuMP.Model = SOSModel())


"""
function chordal_sos(p::APL; model::JuMP.Model = SOSModel())
    H, cliques = chordal_csp_graph(p)
    degree_p = MP.maxdegree(p)
    for clique in cliques
        vars = Tuple(unique!(sort!(clique, rev = true)))
        mvec = MP.monomials(vars, 0:div(degree_p,2))
        pp = @variable model variable_type=SOSPoly(mvec)
        p = p - pp
    end
    @constraint model p == 0
    return model
end

"""
    chordal_putinar(
                         p::APL, 
                         degree::Int,
                         K::AbstractBasicSemialgebraicSet;
                         model::JuMP.Model = SOSModel()
                        )


"""
function chordal_putinar(
                         p::APL, 
                         degree::Int,
                         K::AbstractBasicSemialgebraicSet;
                         model::JuMP.Model = SOSModel()
                        )

    if K isa FullSpace
        println("Unbounded domain. Ignoring degree = $degree .")
        model = chordal_sos(p, model)
    else
        H, cliques = chordal_csp_graph(p, K)

        for clique in cliques

            vars = Tuple(unique!(sort!(clique, rev = true)))
            mvec = MP.monomials(vars, 0:MP.maxdegree(p))
            pp = @variable model variable_type=Poly(mvec) 

            p = p - pp


            Ki = BasicSemialgebraicSet{Float64, typeof(p)}()
            if K isa BasicSemiAlgebraicSet
                for ineq in inequalities(K)
                    if CEG.contains(clique, variables(ineq))
                        addinequality!(Ki,ineq)
                    end
                end
            end

            for eq in equalities(K)
                if CEG.contains(clique, variables(eq))
                    deg = degree - MP.maxdegree(eq)
                    if deg >= 0
                        mvec = MP.monomials(vars, 0:deg)
                        ppi = @variable model variable_type=Poly(mvec)
                        pp = pp - ppi*eq
                    end
                end
            end

            @constraint model pp in SOSCone() domain = Ki maxdegree = degree
        end

        @constraint model p == 0
    end
    return model
end

