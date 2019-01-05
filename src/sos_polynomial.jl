export DSOSCone, SDSOSCone, SOSCone

export CoDSOSCone, CoSDSOSCone, CoSOSCone
export SOSMatrixCone
export getslack, certificate_monomials, lagrangian_multipliers

struct DSOSCone <: PolyJuMP.PolynomialSet end
struct CoDSOSCone <: PolyJuMP.PolynomialSet end
_varconetype(::DSOSCone) = DSOSPoly
_nococone(::CoDSOSCone) = DSOSCone()

struct SDSOSCone <: PolyJuMP.PolynomialSet end
struct CoSDSOSCone <: PolyJuMP.PolynomialSet end
_varconetype(::SDSOSCone) = SDSOSPoly
_nococone(::CoSDSOSCone) = SDSOSCone()

struct SOSCone <: PolyJuMP.PolynomialSet end
struct CoSOSCone <: PolyJuMP.PolynomialSet end
_varconetype(::SOSCone) = SOSPoly
_nococone(::CoSOSCone) = SOSCone()

const SOSLikeCones = Union{DSOSCone, SDSOSCone, SOSCone}
const CoSOSLikeCones = Union{CoDSOSCone, CoSDSOSCone, CoSOSCone}
const SOSSubCones = Union{CoSOSLikeCones, SOSLikeCones}

struct SOSPolynomialSet{DT <: AbstractSemialgebraicSet,
                        CT <: SOSSubCones,
                        BT <: PolyJuMP.AbstractPolynomialBasis,
                        MT <: AbstractMonomial,
                        MVT <: AbstractVector{MT},
                        NPT <: Tuple} <: MOI.AbstractVectorSet
    domain::DT
    cone::CT
    basis::Type{BT}
    monomials::MVT
    newton_polytope::NPT
    mindegree::Int
    maxdegree::Int
end
