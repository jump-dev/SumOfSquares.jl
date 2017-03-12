export createpoly, createnonnegativepoly

function createpoly{C}(m::JuMP.Model, monotype::Symbol, x::MonomialVector{C})
    if monotype == :Default
        monotype = :Classic
    end
    gram = monotype == :Gram
    if gram
        Z = (sum(x)^2).x
    else
        Z = x
    end
    Polynomial{C, JuMP.Variable}((i) -> Variable(m, -Inf, Inf, :Cont), Z)
end
createpoly(m::JuMP.Model, monotype::Symbol, x::Vector) = createpoly(m, monotype, MonomialVector(x))

function createnonnegativepoly{C}(m::JuMP.Model, monotype::Symbol, x::MonomialVector{C})
    if isempty(x)
        # Need MultivariatePolynomials v0.0.2
        #zero(MatPolynomial{C, JuMP.Variable})
        MatPolynomial(JuMP.Variable[], x)
    else
        if monotype == :Default
            monotype = :Gram
        end
        gram = monotype == :Gram
        if gram
            Z = x
        else
            Z = getmonomialsforcertificate(x)
        end
        p = MatPolynomial{C, JuMP.Variable}((i, j) -> Variable(m, -Inf, Inf, :Cont), Z)
        push!(m.varCones, (:SDP, p.Q[1].col:p.Q[end].col))
        if !gram
            # The coefficients of a monomial not in Z do not all have to be zero, only their sum
            addpolyeqzeroconstraint(m, removemonomials(Polynomial(p), x))
        end
        p
    end
end
createnonnegativepoly(m::JuMP.Model, monotype::Symbol, x::Vector) = createnonnegativepoly(m, monotype, MonomialVector(x))
