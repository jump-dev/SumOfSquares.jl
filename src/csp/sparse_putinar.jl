using MultivariatePolynomials
const APL = AbstractPolynomialLike

"""
    sparse_putinar(objective, equalities, inequalities, degree)

"""
function sparse_putinar(objective::T, equalities::Vector{T}, inequalities::Vector{T}, degree::Int) where T <: APL
    H, cliques = chordal_csp_graph(objective, [equalities, inequalities])
