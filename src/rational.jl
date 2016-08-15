import Base./

immutable Rational{S,T}
  num::TermContainer{S}
  den::TermContainer{T}
end

function (/){S,T}(num::TermContainer{S}, den::TermContainer{T})
  Rational{S,T}(num, den)
end
