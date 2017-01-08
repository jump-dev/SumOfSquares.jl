# Adapted from:
# SOSDEMO5 --- Upper bound for the structured singular value mu
# Section 3.5 of SOSTOOLS User's Manual

facts("SOSDEMO5") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
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
        A = Vector{VecPolynomial{Float64}}(4)

        for (gam, expected) in [(0.8723, :Infeasible), (0.8724, :Optimal)]
            Z = monomials(x, 1)
            for i = 1:4
                H = M[i,:]*M[i,:]' - (gam^2)*sparse([i],[i],[1],4,4)
                H = [real(H) -imag(H); imag(H) real(H)]
                A[i] = dot(Z, H*Z)
            end

            m = Model(solver = solver)

            # -- Q(x)'s -- : sums of squares
            # Monomial vector: [x1; ... x8]
            Q = Vector{MatPolynomial{JuMP.Variable}}(4)
            for i = 1:4
                @polyvariable m tmp >= 0 Z
                Q[i] = tmp
            end

            # -- r's -- : constant sum of squares
            Z = monomials(x, 0)
            #r = Matrix{MatPolynomial{JuMP.Variable}}(4,4) # FIXME doesn't work with 1x1 SDP matrix :(
            r = Matrix{JuMP.Variable}(4,4)
            for i = 1:4
                for j = (i+1):4
                    #@polyvariable m tmp >= 0 Z
                    @variable m tmp >= 0
                    r[i,j] = tmp
                end
            end

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

            @polyconstraint m expr >= 0

            status = solve(m)

            # Program is feasible, thus 0.8724 is an upper bound for mu.
            @fact status --> expected
        end
    end
end; end; end
