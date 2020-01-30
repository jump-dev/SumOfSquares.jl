
function monomial_sparsity(f::MP.AbstractPolynomialLike, K::BasicSemialgebraicSet, d::Int)
    @assert isempty(equalities(K))

    vars = unique!(sort!([variables(K)...,variables(f)...], rev=true))
    cons = [1, inequalities(K)...]

    ds = Dict(con =>div(d- maxdegree(con)+ 1, 2) for con in cons)

    M = [mon for mon in monomials(f)]
    for con in cons
        for mon in monomials(con)
            push!(M, mon)
        end
    end

    unique!(sort!(M, rev=true))
    multiplier_bases = Dict{eltype(cons), Any}()
    G = Dict{eltype(cons), CEG.LabelledGraph{eltype(M)}}(con => CEG.LabelledGraph{eltype(M)}() for con in cons )
    for con in cons
        CEG.add_node!.(G[con], [mon for mon in monomials(vars, 0:ds[con])])
    end
    finish = false
    while !finish
        finish = true
        for con in cons
            for i in 1:CEG.num_nodes(G[con].graph)
                for j in i+1:CEG.num_nodes(G[con].graph)
                    if !(j in CEG.neighbors(G[con].graph, i))&&!isempty(intersect([G[con].int2n[i]*G[con].int2n[j]*mon for mon in monomials(con)], M))
                        finish = false
                        CEG.add_edge!(G[con].graph, i, j)
                        for mon in monomials(con)
                            push!(M, mon*G[con].int2n[i]*G[con].int2n[j])
                        end
                    end    
                end
            end
            if !finish
                G[con], multiplier_bases[con] = CEG.chordal_extension(G[con], CEG.GreedyFillIn())
            end
        end
    end
    return multiplier_bases
end

function mono_sparse_certificate(m::Model, f::MP.AbstractPolynomialLike, K::BasicSemialgebraicSet, d::Int)
    multipliers = monomial_sparsity(f, K, d)
    p = f
    for (con, mvs) in multipliers
        for mv in mvs
            mult = @variable(m, [1], SOSPoly(mv)) 
            p -= con*mult[1]
        end
    end
    @constraint(m, p == 0)
    return m
end
