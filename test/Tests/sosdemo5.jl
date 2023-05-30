# Adapted from:
# SOSDEMO5 --- Upper bound for the structured singular value mu
# Section 3.5 of SOSTOOLS User's Manual

using SparseArrays

function sosdemo5_test(optimizer, config::MOI.Test.Config, feasible::Bool, γ)
    @polyvar x[1:8]

    # The matrix under consideration
    α = 3 + sqrt(3)
    β = sqrt(3) - 1
    a = sqrt(2 / α)
    b = 1 / sqrt(α)
    c = b
    d = -sqrt(β / α)
    f = (1 + im) * sqrt(1 / (α * β))
    U = [a 0; b b; c im*c; d f]
    V = [0 a; b -b; c -im*c; -im*f -d]
    M = U * V'

    Z = monomials(x, 1)

    function build_A(i)
        H = M[i, :] * M[i, :]' - (γ^2) * sparse([i], [i], [1], 4, 4)
        H = [real(H) -imag(H); imag(H) real(H)]
        return dot(Z, H * Z)
    end
    A = build_A.(1:4)

    model = _model(optimizer)
    PolyJuMP.setpolymodule!(model, SumOfSquares)

    # -- Q(x)'s -- : sums of squares
    # Monomial vector: [x1; ... x8]
    Q = Vector{
        GramMatrix{
            JuMP.VariableRef,
            monomial_type(x[1]),
            monomial_vector_type(x[1]),
        },
    }(
        undef,
        4,
    )
    @variable model Q[1:4] SOSPoly(Z)

    # -- r's -- : constant sum of squares
    Z = monomials(x, 0)
    #r = Matrix{GramMatrix{JuMP.VariableRef}}(4,4) # FIXME doesn't work with 1x1 SDP matrix :(
    @variable model r[i = 1:4, j = (i+1):4] >= 0

    # Constraint : -sum(Qi(x)*Ai(x)) - sum(rij*Ai(x)*Aj(x)) + I(x) >= 0
    expr = 0
    # Adding term
    for i in 1:4
        expr -= A[i] * Q[i]
    end
    for i in 1:4
        for j in (i+1):4
            expr -= A[i] * A[j] * r[i, j]
        end
    end
    # Constant term: I(x) = -(x1^4 + ... + x8^4)
    I = -sum(x .^ 4)
    expr = expr + I

    @constraint model expr >= 0

    JuMP.optimize!(model)

    # Program is feasible, thus 0.8724 is an upper bound for mu.
    if feasible
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    else
        @test JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
    end
end
const γ_opt = 0.8723
function sosdemo5_infeasible_test(optimizer, config, ε = 1e-2)
    return sosdemo5_test(optimizer, config, false, γ_opt - ε)
end
function sosdemo5_feasible_test(optimizer, config, ε = 1e-2)
    return sosdemo5_test(optimizer, config, true, γ_opt + ε)
end
sd_tests["sosdemo5_infeasible"] = sosdemo5_infeasible_test
sd_tests["sosdemo5_feasible"] = sosdemo5_feasible_test
