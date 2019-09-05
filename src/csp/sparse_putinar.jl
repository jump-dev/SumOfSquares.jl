export chordal_sos, chordal_putinar

"""
    chordal_sos(p::APL; model::JuMP.Model = SOSModel())


"""
function chordal_sos(p::APL; model::JuMP.Model = SOSModel())
    H, cliques = chordal_csp_graph(p)
    degree_p = MP.maxdegree(p)
    for clique in cliques
        vars = Tuple(unique!(sort!(clique, rev = true)))
        mvec = MP.monomials(vars, 0:degree_p)
        pp = @variable model [1] Poly(mvec)
        p = p - pp[1]
        @constraint model pp[1] in SOSCone()
    end
    @constraint model p == 0
    return model
end

"""
chordal_putinar(p::T, degree::Int; equalities = T[], inequalities = T[],  model = SOSModel)


"""
function chordal_putinar(
                         p::T, 
                         degree::Int;
                         equalities = T[], 
                         inequalities = T[], 
                         model::JuMP.Model = SOSModel()
                        ) where T <: APL
    if isempty([equalities, inequalities])
        model = chordal_sos(p, model)
    else
        H, cliques = chordal_csp_graph(p, [equalities..., inequalities...])
        println(cliques)
        for clique in cliques
            vars = Tuple(unique!(sort!(clique, rev = true)))
            mvec = MP.monomials(vars, 0:degree)
            pp = @variable model [1] Poly(mvec) 
            p = p - pp[1]
            for eq in equalities
                if CEG.contains(clique, variables(eq))
                    deg = degree - MP.maxdegree(eq)
                    if deg >= 0
                        vars = Tuple(unique!(sort!(clique, rev = true)))
                        mvec = MP.monomials(vars, 0:deg)
                        ppi = @variable model [1] Poly(mvec)
                        pp = pp - ppi[1]*eq
                    end
                end
            end
            for ineq in inequalities
                if CEG.contains(clique, variables(ineq))
                    println(clique)
                    println(ineq)
                    println()
                    degree_s = 2*div( degree - MP.maxdegree(ineq) + 1, 2)
                    if degree_s >= 0
                        vars = Tuple(unique!(sort!(clique, rev = true)))
                        mvec = MP.monomials(vars, 0:degree_s)
                        if degree_s == 0
                            ppi = @variable model [1] 
                            @constraint model ppi[1]>=0
                        else
                            ppi = @variable model [1] Poly(mvec)
                            @constraint model ppi[1] in SOSCone()
                        end
                        pp = pp - ppi[1]*ineq
                    end
                end
            end
        end
        @constraint model p == 0
    end
    return model
end
