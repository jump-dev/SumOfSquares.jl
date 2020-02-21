"""
Classical test functions for sparse polynomial optimization. See e.g. [WSMM06].

[WSMM06] Waki, Hayato, Sunyoung Kim, Masakazu Kojima, and Masakazu Muramatsu.
Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity.
SIAM Journal on Optimization 17, no. 1 (2006): 218-242.
"""
function chained_singular(n::Int)
    # clique size 3
    @assert mod(n, 4)==0
    @polyvar x[1:n]
    return sum([2*i-1 for i in 1:Int(n/2-1)]) do j
        (x[j] + 10*x[j+1])^2 + 5*(x[j+2] - x[j+3])^2 + (x[j+1] - 2*x[j+2])^4 + 10*(x[j] - x[j+3])^4
    end
end

function broyden_banded(n::Int)
    # clique size 7
    @polyvar x[1:n]
    return sum(1:n) do i
        ( x[i]*(2+5*x[i]^2) + 1 + (1 + x[i])*x[i] - sum( (1+x[j])*x[j] for j = maximum([1, i-5]):minimum([n, i+1]) ) )^2
    end
end

function broyden_tridiagonal(n::Int)
    # clique size 3
    @polyvar x[1:n]
    return (( 3 - 2*x[1])*x[1] -2*x[2] +1 )^2 + sum( ((3 - 2*x[i])*x[i] - x[i-1] - 2*x[i+1] + 1)^2 for i = 2:n-1)
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

function sos_lower_bound(p, factory, sparsity::Sparsity)
    model = Model(factory)
    @variable(model, t)
    @objective(model, Max, t)
    @constraint(model, p - t in SOSCone(), sparsity=sparsity)
    optimize!(model)
    println(termination_status(model))
    println(objective_value(model))
    return model
end
