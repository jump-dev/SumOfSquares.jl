export polytype, createpoly, SOSPoly

function JuMP.getvalue(p::MatPolynomial{JuMP.Variable})
    MatPolynomial(map(getvalue, p.Q), p.x)
end

struct SOSPoly{MT, MV} <: PolyJuMP.AbstractPoly
    x::MV
end
SOSPoly{MT}(x::MV) where {MT, MV} = SOSPoly{MT, MV}(x)
SOSPoly(x::MV) where MV = SOSPoly{:Default}(x)

const PosPoly{MT, MV} = Union{SOSPoly{MT, MV}, Poly{true, MT, MV}}

polytype(m::JuMP.Model, p) = _polytype(m, p, p.x)
polytype(m::JuMP.Model, p, X::AbstractVector) = _polytype(m, p, monovec(X))

# Free polynomial

_polytype(m::JuMP.Model, ::Poly{false}, x::AbstractVector{MT}) where MT<:AbstractMonomial = polynomialtype(MT, JuMP.Variable)

# x should be sorted and without duplicates
function _createpoly(m::JuMP.Model, ::Poly{false}, x::AbstractVector{<:AbstractMonomial}, category::Symbol)
    polynomial((i) -> Variable(m, -Inf, Inf, category), x)
end
function createpoly(m::JuMP.Model, p::Union{Poly{false, :Default}, Poly{false, :Classic}}, category::Symbol)
    _createpoly(m, p, monovec(p.x), category)
end
function createpoly(m::JuMP.Model, p::Poly{false, :Gram}, category::Symbol)
    _createpoly(m, p, monomials(sum(p.x)^2), category)
end

# Sum-of-Squares polynomial

_polytype(m::JuMP.Model, ::PosPoly, x::MVT) where {MT<:AbstractMonomial, MVT<:AbstractVector{MT}} = MatPolynomial{JuMP.Variable, MT, MVT}

function _createpoly(m::JuMP.Model, pp::PosPoly, x::AbstractVector{<:AbstractMonomial}, category::Symbol)
    if isempty(x)
        zero(_polytype(m, pp, x))
    else
        p = MatPolynomial{JuMP.Variable}((i, j) -> Variable(m, -Inf, Inf, category), x)
        push!(m.varCones, (:SDP, p.Q[1].col:p.Q[end].col))
        p
    end
end
function createpoly(m::JuMP.Model, p::Union{PosPoly{:Default}, PosPoly{:Gram}}, category::Symbol)
    _createpoly(m, p, monovec(p.x), category)
end
function createpoly(m::JuMP.Model, pp::PosPoly{:Classic}, category::Symbol)
    p = _createpoly(m, pp, getmonomialsforcertificate(pp.x), category)
    # The coefficients of a monomial not in Z do not all have to be zero, only their sum
    addpolyeqzeroconstraint(m, removemonomials(polynomial(p), pp.x))
    p
end
