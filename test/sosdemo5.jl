# Adapted from:
# SOSDEMO5 --- Upper bound for the structured singular value mu
# Section 3.5 of SOSTOOLS User's Manual

@testset "SOSDEMO5 with $(factory.constructor)" for factory in sdp_factories
    isscs(factory) && continue
    iscsdp(factory) && continue # See https://github.com/JuliaOpt/SumOfSquares.jl/issues/52
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

    @testset "with γ=$γ it should be $(feasible ? "feasible" : "infeasible")" for (γ, feasible) in ((0.8722, false), (0.8724, true))
        Z = monomials(x, 1)

        function build_A(i)
            H = M[i,:]*M[i,:]' - (γ^2)*sparse([i],[i],[1],4,4)
            H = [real(H) -imag(H); imag(H) real(H)]
            dot(Z, H*Z)
        end
        A = build_A.(1:4)

        m = SOSModel(factory)

        # -- Q(x)'s -- : sums of squares
        # Monomial vector: [x1; ... x8]
        Q = Vector{GramMatrix{JuMP.VariableRef, monomialtype(x[1]),
                                 monovectype(x[1])}}(undef, 4)
        @variable m Q[1:4] SOSPoly(Z)

        # -- r's -- : constant sum of squares
        Z = monomials(x, 0)
        #r = Matrix{GramMatrix{JuMP.VariableRef}}(4,4) # FIXME doesn't work with 1x1 SDP matrix :(
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

        JuMP.optimize!(m)

        # Program is feasible, thus 0.8724 is an upper bound for mu.
        if feasible
            @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        else
            @test JuMP.dual_status(m) == MOI.INFEASIBILITY_CERTIFICATE
        end
    end
end
