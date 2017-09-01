# Adapted from:
# SOSDEMO5 --- Upper bound for the structured singular value mu
# Section 3.5 of SOSTOOLS User's Manual

@testset "SOSDEMO5 with $solver" for solver in sdp_solvers
    if !isscs(solver)
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

        # Constructing A(x)'s
        A = Vector{Polynomial{true, Float64}}(4)

        for (gam, expected) in [(0.8723, :Infeasible), (0.8724, :Optimal)]
            Z = monomials(x, 1)
            for i = 1:4
                H = M[i,:]*M[i,:]' - (gam^2)*sparse([i],[i],[1],4,4)
                H = [real(H) -imag(H); imag(H) real(H)]
                A[i] = dot(Z, H*Z)
            end

            m = SOSModel(solver = solver)

            # -- Q(x)'s -- : sums of squares
            # Monomial vector: [x1; ... x8]
            Q = Vector{MatPolynomial{JuMP.Variable, monomialtype(x[1]), monovectype(x[1])}}(4)
            @variable m Q[1:4] >= 0 Poly{true}(Z)

            # -- r's -- : constant sum of squares
            Z = monomials(x, 0)
            #r = Matrix{MatPolynomial{JuMP.Variable}}(4,4) # FIXME doesn't work with 1x1 SDP matrix :(
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

            status = solve(m)

            # Program is feasible, thus 0.8724 is an upper bound for mu.
            @test status == expected
        end
    end
end
