function decompose(p::MP.AbstractPolynomial, tol)
    vars = MP.effective_variables(p)
    if length(vars) != 1
        error("$p is not univariate")
    end
    x = first(vars)
    lead = MP.leadingcoefficient(p)
    if !isone(lead)
        p = p / lead
    end
    deg = MP.maxdegree(p)
    if isodd(deg)
        return
    end
    d = div(deg, 2)
    companion = zeros(2d, 2d)
    for i in 0:(2d-1)
        if i > 0
            companion[i + 1, i] = 1.0
        end
        companion[i + 1, end] = -MP.coefficient(p, x^i)
    end
    #display(companion)
    F = LinearAlgebra.schur(complex(companion))
    #display(F.Z * F.T * F.Z')
    σ = vcat(1:2:(2d - 1), 2:2:2d)
    S = LinearAlgebra.ordschur(F, isodd.(1:2d))
    #U = S.Z[σ, σ]
    #S.Z * S.T * S.Z'
    #display(σ)
    #display(S.T)
    #S.T[σ, σ]
    vt = S.Z[1, 1:d]
    I = 1:d
    J = d .+ I
    U11 = S.Z[I, I]
    U22 = S.Z[J, J]
    U12 = S.Z[I, J]
    display(S.Z)
    return U12
    U21 = S.Z[J, I]
    vt = U22[:, 1]
    U = [U11 U12; U21 U22]
    display(S.T)
    display(U * S.T * U')
    display(S.Z * S.T * S.Z')
    display(U * U')
    @show lead
    q = lead * U12' \ vt
    q = vec(vt' * inv(U12))
    @show q
    q1 = lead * x^d + sum([-real(q[i + 1]) * x^i for i in 0:(d-1)])
    q2 = sum([-imag(q[i + 1]) * x^i for i in 0:(d-1)])
    q_1 = Polynomial([-real.(q); 1.0])
    q_2 = Polynomial([-imag.(q); 0.0])
    @show q1^2 + q2^2
    @show q_1^2 + q_2^2
    @show Polynomials.roots(q_1 + im * q_2)
    @show Polynomials.roots(q_1 - im * q_2)
    return q1, q2
end
