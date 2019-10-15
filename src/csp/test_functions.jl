"""
    Classical test functions for sparse polynomial optimization. See e.g. HAYATO  WAKI,  SUNYOUNG  KIM,  MASAKAZU  KOJIMA,AND MASAKAZU MURAMATSU: SUMS OF SQUARES AND SEMIDEFINITE PROGRAMRELAXATIONS FOR POLYNOMIAL OPTIMIZATION PROBLEMSWITH STRUCTURED SPARSITY, SIAM J. OPTIM. Vol. 17, No. 1, pp. 218–242  SIAM J. OPTIM. Vol. 17, No. 1, pp. 218–242 
"""
function chained_singular(n::Int)
    # clique size 3
    @assert mod(n, 4)==0
    @polyvar x[1:n]
    p = sum( (x[j] + 10*x[j+1])^2 + 5*(x[j+2] - x[j+3])^2 + (x[j+1] - 2*x[j+2])^4 + 10*(x[j] - x[j+3])^4 for j in [2*i-1 for i in 1:Int(n/2-1)]) 
    return p
end

function broyden_banded(n::Int)
    # clique size 7
    @polyvar x[1:n]
    p = sum( ( x[i]*(2+5*x[i]^2) + 1 + (1 + x[i])*x[i] - sum( (1+x[j])*x[j] for j = maximum([1, i-5]):minimum([n, i+1]) ) )^2 for i=1:n)
    return p
end

function broyden_tridiagonal(n::Int)
    # clique size 3
    @polyvar x[1:n]
    p = (( 3 - 2*x[1])*x[1] -2*x[2] +1 )^2 + sum( ((3 - 2*x[i])*x[i] - x[i-1] - 2*x[i+1] + 1)^2 for i = 2:n-1)
end

function chained_wood(n::Int)
    # clique size 2
    @assert mod(n, 4)==0
    @polyvar x[1:n]
    p = 1 + sum(100*(x[j+1] - x[j]^2)^2 + (1 - x[j])^2 + 90*(x[j+3] - x[j+2]^2)^2  + (1 - x[j+2])^2 
            + 10*(x[j+1] + x[j+3] - 2)^2 + 0.1*(x[j+1] - x[j+3])^4 for j in [2*i-1 for i in 1:Int(n/2-1)])
    return p
end

function generalized_rosenbrock(n::Int)
    # clique size 2
    @polyvar x[1:n]
    p = 1 + sum( 100*(x[i]-x[i-1]^2)^2 + (1-x[i])^2 for i=2:n)
    return p
end

function dense(p)
    dense_model = SOSModel(with_optimizer(Mosek.Optimizer))
    @variable dense_model t
    @objective dense_model Max t
    @constraint dense_model p-t in SOSCone()
    optimize!(dense_model)
    println(termination_status(dense_model))
    println(objective_value(dense_model))
    return dense_model
end

function sparse(p)
    sparse_model = SOSModel(with_optimizer(Mosek.Optimizer))
    @variable sparse_model t
    @objective sparse_model Max t
    chordal_sos(p-t, model=sparse_model)
    optimize!(sparse_model)
    println(termination_status(sparse_model))
    println(objective_value(sparse_model))
    return sparse_model
end
