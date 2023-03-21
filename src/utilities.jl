# The second type is ignored. For `T = VariableRef`, it should be `Float64` anyway
# since JuMP only supports `Float64`.
_promote_sum(T::Type, ::Type = Float64) = MA.promote_operation(+, T, T)
# `+` is not defined between `MOI.VariableIndex`.
function _promote_sum(::Type{MOI.VariableIndex}, T::Type = Float64)
    return MOI.ScalarAffineFunction{T}
end

_promote_add_mul(T::Type) = MA.promote_operation(MA.add_mul, T, T, T)
function _promote_add_mul(::Type{MOI.VariableIndex})
    return MOI.ScalarQuadraticFunction{Float64}
end
function MP.polynomial(p::GramMatrix{MOI.VariableIndex,B,U}) where {B,U}
    Q = convert(Vector{U}, p.Q.Q)
    return MP.polynomial(GramMatrix(SymMatrix(Q, p.Q.n), p.basis))
end
#function MP.polynomial(p::GramMatrix{F}) where {F <: MOI.AbstractFunction}
#    MP.polynomial(p, MOI.Utilities.promote_operation(+, Float64, F, F))
#end

function primal_value(model, p::GramMatrix{MOI.VariableIndex})
    Q = MOI.get(model, MOI.VariablePrimal(), p.Q.Q)
    return GramMatrix(SymMatrix(Q, p.Q.n), p.basis)
end

_complex(::Type{T}, ::Type) where {T} = T
function _complex(
    ::Type{T},
    ::Type{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle},
) where {T}
    return Complex{T}
end
