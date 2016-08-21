import Base.+, Base.-, Base.*, Base./

immutable Rational{S,T}
  num::TermContainer{S}
  den::TermContainer{T}
end

function (/){S,T}(num::TermContainer{S}, den::TermContainer{T})
  Rational{S,T}(num, den)
end

function (+)(r::Rational, s::Rational)
  (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function (+)(p::VecPolynomial, r::Rational)
  (p*r.den + r.num) / r.den
end
(+)(r::Rational, p::VecPolynomial) = p + r
function (-)(r::Rational, s::Rational)
  (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
(-)(p::PolyType, s::Rational) = (p * s.den - s.num) / s.den
(-)(s::Rational, p::PolyType) = (s.num - p * s.den) / s.den

function (*)(p::TermContainer, r::Rational)
  if p == r.den
    r.num
  else
    (p * r.num) / r.den
  end
end
(*)(r::Rational, p::Term) = p * r
(*)(r::Rational, p::VecPolynomial) = p * r

zero(r::Rational) = zero(r.num)
zero{T<:Rational}(::Type{T}) = zero(VecPolynomial)
zero{T}(::Type{Rational{T}}) = zero(VecPolynomial{T})
