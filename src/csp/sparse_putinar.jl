export sos, chordal_sos, putinar, chordal_putinar


"""
    sos(p::APL; model::JuMP.Model = SOSModel())


"""
function sos(p::APL; model::JuMP.Model = SOSModel())
    degree_p = MP.maxdegree(p)
    vars = Tuple(unique!(sort!(variables(p), rev = true)))
    mvec = MP.monomials(vars, 0:div(degree_p,2))
    pp = @variable model variable_type=SOSPoly(mvec)
    @constraint model p-pp == 0
    return model
end

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
    putinar(
                         p::APL, 
                         degree::Int,
                         K::AbstractBasicSemialgebraicSet;
                         model::JuMP.Model = SOSModel()
                        )


"""
function putinar(
                         p::APL, 
                         degree::Int,
                         K::AbstractBasicSemialgebraicSet;
                         model::JuMP.Model = SOSModel()
                        )

    if K isa FullSpace
        println("Unbounded domain. Ignoring degree = $degree .")
        model = sos(p, model)
    else
        vars = variables(p)
        if K isa BasicSemialgebraicSet
            for ineq in inequalities(K)
                append!(vars, variables(ineq))
            end
        end
        for eq in equalities(K)
            append!(vars, variables(eq))
        end
        vars = Tuple(unique!(sort!(vars, rev = true)))

        while !isempty(equalities(K))
            eq = pop!(equalities(K))
            deg = degree - MP.maxdegree(eq)
            if deg >= 0
                 mvec = MP.monomials(vars, 0:deg)
                 ppi = @variable model variable_type=Poly(mvec)
                 println(ppi)
                 p = p - ppi*eq
            end
        end
        @constraint model p in SOSCone() domain = K maxdegree = degree
    end
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
            mvec = MP.monomials(vars, 0:degree)
            pp = @variable model variable_type=Poly(mvec) 

            p = p - pp

            dummy = sum(rand()*mvec[i] for i=1:length(mvec))
            Ki = BasicSemialgebraicSet{Float64, typeof(dummy)}()

            if K isa BasicSemialgebraicSet
                for ineq in inequalities(K)
                    if CEG.contains(clique, variables(ineq))
                        addinequality!(Ki, ineq)
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

