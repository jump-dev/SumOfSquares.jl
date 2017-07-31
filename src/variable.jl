export polytype, createpoly, createnonnegativepoly

polytype{C}(m::JuMP.Model, x::MonomialVector{C}) = Polynomial{C, JuMP.Variable}
polytype(m::JuMP.Model, x::Vector) = polytype(m, MonomialVector(x))
polytype(m::JuMP.Model, p::Poly) = polytype(m, p.x)

function createpoly{C}(m::JuMP.Model, x::MonomialVector{C}, category::Symbol)
    Polynomial{C, JuMP.Variable}((i) -> Variable(m, -Inf, Inf, category), x)
end
createpoly(m::JuMP.Model, x::Vector, category::Symbol) = createpoly(m, MonomialVector(x), category)
function createpoly(m::JuMP.Model, p::Union{Poly{false, :Default}, Poly{false, :Classic}}, category::Symbol)
    createpoly(m, p.x, category)
end
function createpoly(m::JuMP.Model, p::Poly{false, :Gram}, category::Symbol)
    createpoly(m, (monomials(sum(p.x)^2)), category)
end

nonnegativepolytype{C}(m::JuMP.Model, x::MonomialVector{C}) = MatPolynomial{C, JuMP.Variable}
nonnegativepolytype(m::JuMP.Model, x::Vector) = nonnegativepolytype(m, MonomialVector(x))
nonnegativepolytype(m::JuMP.Model, p::Poly) = nonnegativepolytype(m, p.x)

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
function createnonnegativepoly(m::JuMP.Model, p::Union{Poly{true, :Default}, Poly{true, :Gram}}, category::Symbol)
    createnonnegativepoly(m, p.x, category)
end
function createnonnegativepoly(m::JuMP.Model, p::Poly{true, :Classic}, category::Symbol)
    p = createnonnegativepoly(m, getmonomialsforcertificate(p), category)
    # The coefficients of a monomial not in Z do not all have to be zero, only their sum
    addpolyeqzeroconstraint(m, removemonomials(Polynomial(p), monomials(p)), FullSpace())
    p
end
