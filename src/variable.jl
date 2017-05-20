export createpoly, createnonnegativepoly

function createpoly{C}(m::JuMP.Model, x::MonomialVector{C}, category::Symbol)
    Polynomial{C, JuMP.Variable}((i) -> Variable(m, -Inf, Inf, category), x)
end
createpoly(m::JuMP.Model, x::Vector, category::Symbol) = createpoly(m, MonomialVector(x), category)
function createpoly(m::JuMP.Model, p::Union{Poly{:Default}, Poly{:Classic}}, category::Symbol)
    createpoly(m, p.x, category)
end
function createpoly(m::JuMP.Model, p::Union{:Gram}, category::Symbol)
    createpoly(m, (sum(p.x)^2).x, category)
end
function createpoly{MT}(m::JuMP.Model, p::Poly{MT}monotype::Symbol, (sum(p.x)^2).x, category::Symbol)
    if monotype == :Default
        monotype = :Classic
    end
    gram = monotype == :Gram
    if gram
        Z = (sum(x)^2).x
    else
        Z = x
    end
    Polynomial{C, JuMP.Variable}((i) -> Variable(m, -Inf, Inf, category), Z)
end

function createnonnegativepoly{C}(m::JuMP.Model, x::MonomialVector{C}, category::Symbol)
    if isempty(x)
        # Need MultivariatePolynomials v0.0.2
        #zero(MatPolynomial{C, JuMP.Variable})
        MatPolynomial(JuMP.Variable[], x)
    else
        p = MatPolynomial{C, JuMP.Variable}((i, j) -> Variable(m, -Inf, Inf, category), x)
        push!(m.varCones, (:SDP, p.Q[1].col:p.Q[end].col))
        p
    end
end
function createnonnegativepoly(m::JuMP.Model, x::Vector, category::Symbol)
    createnonnegativepoly(m, MonomialVector(x), category)
end
function createnonnegativepoly(m::JuMP.Model, p::Union{Poly{:Default}, Poly{:Gram}}, category::Symbol)
    createnonnegativepoly(m, p.x, category)
end
function createnonnegativepoly(m::JuMP.Model, p::Poly{:Classic}, category::Symbol)
    p = createnonnegativepoly(m, getmonomialsforcertificate(p.x), category)
    # The coefficients of a monomial not in Z do not all have to be zero, only their sum
    addpolyeqzeroconstraint(m, removemonomials(Polynomial(p), p.x))
    p
end
