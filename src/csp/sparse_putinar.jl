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
chordal_putinar(p::T, degree::Int; equalities = T[], inequalities = T[],  model = SOSModel)


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
    elseif K isa AbstractAlgebraicSet
        H, cliques = chordal_csp_graph(p, equalities(K))
        for clique in cliques
            vars = Tuple(unique!(sort!(clique, rev = true)))
            mvec = MP.monomials(vars, 0:MP.maxdegree(p))
            pp = @variable model variable_type=Poly(mvec) 
            p = p - pp
            for eq in equalities
                if CEG.contains(clique, variables(eq))
                    deg = degree - MP.maxdegree(eq)
                    if deg >= 0
                        mvec = MP.monomials(vars, 0:deg)
                        ppi = @variable model variable_type=Poly(mvec)
                        pp = pp - ppi[1]*eq
                    end
                end
            end
            @constraint model pp in SOSCone()
        end
        @constraint model p ==0
    else
        if isempty(equalities(K))
            constraints =  [equalities(K)..., inequalities(K)...]
        else
            constraints = inequalities(K)
        end
        
        H, cliques = chordal_csp_graph(p,constraints)

        for clique in cliques
            vars = Tuple(unique!(sort!(clique, rev = true)))
            mvec = MP.monomials(vars, 0:MP.maxdegree(p))
            pp = @variable model variable_type=Poly(mvec) 
            p = p - pp
            for eq in equalities(K)
                if CEG.contains(clique, variables(eq))
                    deg = degree - MP.maxdegree(eq)
                    if deg >= 0
                        mvec = MP.monomials(vars, 0:deg)
                        ppi = @variable model variable_type=Poly(mvec)
                        pp = pp - ppi[1]*eq
                    end
                end
            end
            for ineq in inequalities(K)
                if CEG.contains(clique, variables(ineq))
                    degree_s = div(degree - MP.maxdegree(ineq) + 1, 2)
                    if degree_s >= 0
                        mvec = MP.monomials(vars, 0:degree_s)
                        ppi = @variable model variable_type=SOSPoly(mvec)   
                        pp = pp - ppi*ineq
                    end
                end
            end
            @constraint model pp in SOSCone()
        end
        @constraint model p == 0
    end
    return model
end

