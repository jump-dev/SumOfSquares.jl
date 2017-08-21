export polytype, createpoly, SOSPoly

function JuMP.getvalue{C}(p::MatPolynomial{JuMP.Variable})
    MatPolynomial(map(getvalue, p.Q), p.x)
end

struct SOSPoly{MT, MV} <: PolyJuMP.AbstractPoly
    x::MV
end
SOSPoly{MT}(x::MV) where {MT, MV} = SOSPoly{MT, MV}(x)
SOSPoly(x::MV) where MV = SOSPoly{:Default}(x)

const PosPoly{MT, MV} = Union{SOSPoly{MT, MV}, Poly{true, MT, MV}}

function _createpoly(m::JuMP.Model, p, x::Vector, category::Symbol)
    _createpoly(m, p, MonomialVector(x), category)
end
_polytype(m::JuMP.Model, p, x::Vector) = _polytype(m, p, MonomialVector(x))
polytype(m::JuMP.Model, p) = _polytype(m, p, p.x)

# Free polynomial

_polytype{C}(m::JuMP.Model, ::Poly{false}, x::MonomialVector{C}) = Polynomial{C, JuMP.Variable}

function _createpoly{C}(m::JuMP.Model, ::Poly{false}, x::MonomialVector{C}, category::Symbol)
    Polynomial{C, JuMP.Variable}((i) -> Variable(m, -Inf, Inf, category), x)
end
function createpoly(m::JuMP.Model, p::Union{Poly{false, :Default}, Poly{false, :Classic}}, category::Symbol)
    _createpoly(m, p, p.x, category)
end
function createpoly(m::JuMP.Model, p::Poly{false, :Gram}, category::Symbol)
    _createpoly(m, p, (sum(p.x)^2).x, category)
end

# Sum-of-Squares polynomial

_polytype{C}(m::JuMP.Model, ::PosPoly, x::MonomialVector{C}) = MatPolynomial{C, JuMP.Variable}

function _createpoly{C}(m::JuMP.Model, ::PosPoly, x::MonomialVector{C}, category::Symbol)
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
function createpoly(m::JuMP.Model, p::Union{PosPoly{:Default}, PosPoly{:Gram}}, category::Symbol)
    _createpoly(m, p, p.x, category)
end
function createpoly(m::JuMP.Model, p::PosPoly{:Classic}, category::Symbol)
    p = _createpoly(m, p, getmonomialsforcertificate(p.x), category)
    # The coefficients of a monomial not in Z do not all have to be zero, only their sum
    addpolyeqzeroconstraint(m, removemonomials(Polynomial(p), p.x))
    p
end
