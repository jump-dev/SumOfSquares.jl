facts("Monomial selection for certificate") do
    @polyvar x y
    @fact_throws ErrorException getmonomialsforcertificate([x*y, y^2], :Sparse)
    @fact getmonomialsforcertificate([x*y, y^2]) --> MonomialVector([y])
end

facts("Random SOS should be SOS") do
    for solver in sdp_solvers
        context("With solver $(typeof(solver))") do
            @polyvar x y
            x = [Monomial(), x, y, x^2, y^2, x*y]
            @fact_throws ArgumentError randsos(x, monotype=:Unknown)
            @fact typeof(x) --> Vector{Monomial}
            for i in 1:10
                for monotype in [:Classic, :Gram]
                    p = randsos(x, monotype=monotype)

                    m = Model(solver = solver)

                    @polyconstraint m p >= 0

                    status = solve(m)

                    @fact status --> :Optimal
                end
            end
end; end; end
