export createpoly, createnonnegativepoly

function createpoly(m::JuMP.Model, ::SOS, monotype::Symbol, x::Union{MonomialVector,Vector})
    if monotype == :Default
        monotype = :Classic
    end
    gram = monotype == :Gram
    if gram
        Z = dot(x, x).x
    else
        Z = x
    end
    VecPolynomial{JuMP.Variable}((i) -> Variable(m, -Inf, Inf, :Cont), Z)
end

function createnonnegativepoly(m::JuMP.Model, ::SOS, monotype::Symbol, x::Union{MonomialVector,Vector})
    if monotype == :Default
        monotype = :Gram
    end
    gram = monotype == :Gram
    if gram
        Z = x
    else
        Z = getmonomialsforcertificate(x)
    end
    p = MatPolynomial{JuMP.Variable}((i,j) -> Variable(m, -Inf, Inf, :Cont), Z)
    push!(m.varCones, (:SDP, p.Q[1].col:p.Q[end].col))
    if !gram
        # The coefficients of a monomial not in Z do not all have to be zero, only their sum
        addpolyeqzeroconstraint(m, removemonomials(VecPolynomial(p), x))
    end
    p
end
