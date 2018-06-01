export DSOSPoly, SDSOSPoly, SOSPoly

function JuMP.getvalue(p::MatPolynomial{JuMP.Variable})
    MatPolynomial(map(getvalue, p.Q), p.x)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{MS, MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: PolyJuMP.AbstractPoly
            x::MV
        end
        $poly{MS}(x::AbstractVector{MT}) where {MS, MT} = $poly{MS, MT, typeof(x)}(x)
        $poly{MS}(x) where MS = $poly{MS}(monovec(x))
        $poly(x) = $poly{:Default}(x)
    end
end

const PosPoly{MT, MV} = Union{DSOSPoly{MT, MV}, SDSOSPoly{MT, MV}, SOSPoly{MT, MV}}

JuMP.variabletype(m::JuMP.Model, p::PosPoly) = PolyJuMP.polytype(m, p, p.x)
PolyJuMP.polytype(m::JuMP.Model, ::PosPoly, x::MVT) where {MT<:AbstractMonomial, MVT<:AbstractVector{MT}} = MatPolynomial{JuMP.Variable, MT, MVT}

function _constraintmatpoly!(m, p, ::Union{SOSPoly, Poly{true}})
    push!(m.varCones, (:SDP, [p.Q[i, j].col for i in 1:size(p.Q, 1) for j in i:size(p.Q, 2)]))
end
function _constraintmatpoly!(m, p, ::DSOSPoly)
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
function _createpoly(m::JuMP.Model, set::PosPoly, x::AbstractVector{<:AbstractMonomial}, category::Symbol)
    p = _matpolynomial(m, x, category)
    if length(x) > 1
        _constraintmatpoly!(m, p, set)
    end
    p
end
function PolyJuMP.createpoly(m::JuMP.Model, p::Union{PosPoly{:Default}, PosPoly{:Gram}}, category::Symbol)
    _createpoly(m, p, monovec(p.x), category)
end
function PolyJuMP.createpoly(m::JuMP.Model, pp::PosPoly{:Classic}, category::Symbol)
    p = _createpoly(m, pp, getmonomialsforcertificate(pp.x), category)
    # The coefficients of a monomial not in Z do not all have to be zero, only their sum
    addpolyconstraint!(m, removemonomials(Polynomial(p), p.x), ZeroPoly(), FullSpace())
    p
end

# Defer other methods to defaults in PolyJuMP
createpoly(args...) = PolyJuMP.createpoly(args...)
