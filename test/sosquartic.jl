 # Test quadratic module when mindegree of polynomial is strictly bigger than 0

@testset "SOSquartic with $(factory.constructor)" for factory in sdp_factories
    @polyvar x  
    m = SOSModel(factory)
    p = x^2*(2-x^2)
    S = @set 1-x^2 >= 0 && x^2 - 0.5 >=0
    soscon = @constraint m p >= 0 domain = S 

    JuMP.optimize!(m)
       
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    Q = gram_matrix(soscon)
    @test mindegree(monomials(Q)) == 0
    @test maxdegree(monomials(Q)) == 4

end
