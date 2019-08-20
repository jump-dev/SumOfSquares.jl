 # Test quadratic module when mindegree of polynomial is strictly bigger than 0

@testset "SOSquartic with $(factory.constructor)" for factory in sdp_factories
    @polyvar x y
    m = SOSModel(factory)
    p = x^4 - y^4 
    S = @set - x^2 - y^4 >=0
    soscon = @constraint m p >= 0 domain = S 

    JuMP.optimize!(m)
       
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    Q = gram_matrix(soscon)
    @test mindegree(monomials(Q)) == 2
    @test maxdegree(monomials(Q)) == 4
end
