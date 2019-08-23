using DynamicPolynomials
#using TypedPolynomials

include("csp_graph.jl")

@polyvar x y z

p = x*y + y*z
q1 = 1 - x^2 - y^2
q2 = 1 - y^2 - z^2

G = csp_graph(variable_type(p), p, [q1, q2])
Matrix(CEG.adjacency_matrix(G))

H, cliques = chordal_csp_graph(variable_type(p), p, [q1, q2])
Matrix(CEG.adjacency_matrix(H))

I, cliques = chordal_csp_graph(p, [q1, q2])
Matrix(CEG.adjacency_matrix(H))


