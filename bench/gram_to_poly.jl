using SumOfSquares
using DynamicPolynomials

function gram_to_poly(n)
    model = MOIU.Model{Float64}()
    @polyvar x
    monos = monomials(x, 0:n)
    q, Q, con_Q = SumOfSquares.add_gram_matrix(model, MOI.PositiveSemidefiniteConeTriangle, monos)
    @time polynomial(q)
    return
end
gram_to_poly(1000)
