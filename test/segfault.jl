using LinearAlgebra
using SparseArrays
using MultivariatePolynomials
using DynamicPolynomials
using JuMP
using SumOfSquares
using CSDP
using Test

solver = CSDP.CSDPSolver(printlevel=0)

    @polyvar x[1:8]

    # The matrix under consideration
    alpha = 3 + sqrt(3);
    beta = sqrt(3) - 1;
    a = sqrt(2/alpha);
    b = 1/sqrt(alpha);
    c = b;
    d = -sqrt(beta/alpha);
    f = (1 + im)*sqrt(1/(alpha*beta));
    U = [a 0; b b; c im*c; d f];
    V = [0 a; b -b; c -im*c; -im*f -d];
    M = U*V';

    for (gam, expected) in [(0.8723, :Infeasible), (0.8724, :Optimal)]
        @show gam
        @show expected
        Z = monomials(x, 1)

        function build_A(i)
            H = M[i,:]*M[i,:]' - (gam^2)*sparse([i],[i],[1],4,4)
            H = [real(H) -imag(H); imag(H) real(H)]
            dot(Z, H*Z)
        end
        A = build_A.(1:4)

        m = SOSModel(solver = solver)

        # -- Q(x)'s -- : sums of squares
        # Monomial vector: [x1; ... x8]
        @variable m Q[1:4] SOSPoly(Z)

        # -- r's -- : constant sum of squares
        Z = monomials(x, 0)
        @variable m r[i=1:4,j=(i+1):4] >= 0

        # Constraint : -sum(Qi(x)*Ai(x)) - sum(rij*Ai(x)*Aj(x)) + I(x) >= 0
        expr = 0
        # Adding term
        for i = 1:4
            expr -= A[i]*Q[i]
        end
        for i = 1:4
            for j = (i+1):4
                expr -= A[i]*A[j]*r[i,j]
            end
        end
        # Constant term: I(x) = -(x1^4 + ... + x8^4)
        I = -sum(x.^4)
        expr = expr + I

        @constraint m expr >= 0

        println("solve")
        status = solve(m)
        println("solved")

        # Program is feasible, thus 0.8724 is an upper bound for mu.
        @test status == expected
    end


