export PseudoExpectation

type PseudoExpectation{T}
  a::Vector{T}
  x::MonomialVector

  function PseudoExpectation(a::Vector{T}, x::MonomialVector)
    if length(a) != length(x)
      error("There should be as many coefficient than monomials")
    end
    zeroidx = Int[]
    for (i,α) in enumerate(a)
      if iszero(α)
        push!(zeroidx, i)
      end
    end
    if !isempty(zeroidx)
      isnz = ones(Bool, length(a))
      isnz[zeroidx] = false
      nzidx = find(isnz)
      a = a[nzidx]
      x = x[nzidx]
    end
    new(a, x)
  end
end

PseudoExpectation{T}(a::Vector{T}, x::MonomialVector) = PseudoExpectation{T}(a, x)

function dot(pexp::PseudoExpectation, p::TermContainer)
  i = 1
  s = 0
  for t in p
    while t.x != pexp.x[i] && i <= length(pexp.x)
      i += 1
    end
    if i > length(pexp.x)
      error("The polynomial $p has a monomial for which the expectation is not known in $pexp")
    end
    s += pexp.a * t.α
    i += 1
  end
  s
end
dot(p::TermContainer, pexp::PseudoExpectation) = dot(pexp, p)
dot(pexp::PseudoExpectation, p::PolyType) = dot(pexp, TermContainer(p))
dot(p::PolyType, pexp::PseudoExpectation) = dot(pexp, TermContainer(p))
(pexp::PseudoExpectation)(p::PolyType) = dot(pexp, p)
