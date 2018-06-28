export DSOSPoly, SDSOSPoly, SOSPoly

function JuMP.getvalue(p::MatPolynomial{JuMP.Variable})
    MatPolynomial(map(getvalue, p.Q), p.x)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{PB<:PolyJuMP.AbstractPolynomialBasis} <: PolyJuMP.AbstractPoly
            polynomial_basis::PB
        end
        $poly(x::AbstractVector{<:MultivariatePolynomials.AbstractPolynomialLike}) = $poly(MonomialBasis(x))
    end
end

const PosPoly{PB} = Union{DSOSPoly{PB}, SDSOSPoly{PB}, SOSPoly{PB}}

JuMP.variabletype(m::JuMP.Model, p::PosPoly) = PolyJuMP.polytype(m, p, p.polynomial_basis)
PolyJuMP.polytype(m::JuMP.Model, ::PosPoly, basis::PolyJuMP.MonomialBasis{MT, MV}) where {MT<:AbstractMonomial, MV<:AbstractVector{MT}} = MatPolynomial{JuMP.Variable, MT, MV}

function _constraintmatpoly!(m, p::MatPolynomial, ::SOSPoly)
    push!(m.varCones, (:SDP, [p.Q[i, j].col for i in 1:size(p.Q, 1) for j in i:size(p.Q, 2)]))
end
function _constraintmatpoly!(m, p::MatPolynomial, ::DSOSPoly)
    n = length(p.x)
    Q = Matrix{JuMP.Variable}(n, n)
    for i in 1:n
        for j in 1:n
            if i == j
                Q[i, j] = p[i, j]
            else
                Q[j, i] = Q[i, j] = Variable(m, -Inf, Inf, :Cont)
                @constraint m Q[i, j] >= p[i, j]
                @constraint m Q[i, j] >= -p[i, j]
            end
        end
    end
    for i in 1:n
        @constraint m 2Q[i, i] >= sum(Q[i, :])
    end
    # If n > 1, this is implied by the constraint but it doesn't hurt to add the variable cone
    # Adding things on varCones makes JuMP think that it is SDP
    # push!(m.varCones, (:NonNeg, map(i -> p[i, i].col, 1:n)))
end
function _matpolynomial(m, x::AbstractVector{<:AbstractMonomial}, category::Symbol)
    if isempty(x)
        zero(JuMP.variabletype(m, SOSPoly(x)))
    else
        if length(x) == 1
            # 1x1 matrix is SDP iff its only entry is nonnegative
            # We handle this case here and do not create any SDP constraint
            lb = 0.
        else
            lb = -Inf
        end
        MatPolynomial{JuMP.Variable}((i, j) -> Variable(m, lb, Inf, category), x)
    end
end
function _createpoly(m::JuMP.Model, set::PosPoly, basis::PolyJuMP.MonomialBasis, category::Symbol)
    p = _matpolynomial(m, basis.monomials, category)
    if length(basis.monomials) > 1
        _constraintmatpoly!(m, p, set)
    end
    p
end
function PolyJuMP.createpoly(m::JuMP.Model, p::PosPoly, category::Symbol)
    _createpoly(m, p, p.polynomial_basis, category)
end

# Defer other methods to defaults in PolyJuMP
createpoly(args...) = PolyJuMP.createpoly(args...)
