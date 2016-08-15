import Base.==, Base.isless, Base.isapprox

# TODO equality should be between name ?
function (==)(x::PolyVar, y::PolyVar)
  x === y
end
isless(x::PolyVar, y::PolyVar) = isless(x.name, y.name)

function isless(x::Vector, y::Vector)
  @assert length(x) == length(y)
  degx = sum(x)
  degy = sum(y)
  if degx != degy
    degx < degy
  else
    for (a, b) in zip(x, y)
      if a < b
        return false
      elseif a > b
        return true
      end
    end
    false
  end
end

# graded lex ordering
function isless(x::Monomial, y::Monomial)
  degx = deg(x)
  degy = deg(y)
  if degx != degy
    degx < degy
  else
    i = j = 1
    while i <= nvars(x)
      @assert j <= nvars(y) # since they have same degree
      if x.vars[i] < y.vars[j]
        if x.z[i] == 0
          i += 1
        else
          return false
        end
      elseif x.vars[i] > y.vars[j]
        if y.z[j] == 0
          j += 1
        else
          return true
        end
      elseif x.z[i] != y.z[j]
        return x.z[i] > y.z[j]
      else
        i += 1
        j += 1
      end
    end
    false # they are equal
  end
end

function isapprox{S,T}(p::VecPolynomial{S}, q::VecPolynomial{T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0)
  i = j = 1
  while i <= length(p.x) || j <= length(q.x)
    lhs, rhs = 0, 0
    if i > length(p.x) || (j <= length(q.x) && q.x[j] < p.x[i])
      @show i
      @show i > length(p.x)
      @show j <= length(q.x) && q.x[j] < p.x[i]
      @show q.x[j]
      @show p.x[i]
      @show q.x[j] < p.x[i]
      rhs = q.a[j]
      j += 1
    elseif j > length(q.x) || p.x[i] < q.x[j]
      @show j
      @show j > length(q.x)
      lhs = p.a[i]
      i += 1
    else
      lhs = p.a[i]
      i += 1
      rhs = q.a[j]
      j += 1
    end
    if !isapprox(lhs, rhs, rtol=rtol, atol=atol)
      @show lhs
      @show rhs
      return false
    end
  end
  true
end

function isapprox{S,T}(p::SOSDecomposition{S}, q::SOSDecomposition{T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0)
  if length(p.ps) != length(q.ps)
    false
  else
    for i in 1:length(p.ps)
      if !isapprox(p.ps[i], q.ps[i], rtol=rtol, atol=atol)
        @show p.ps[i]
        @show q.ps[i]
        return false
      end
    end
    true
  end
end
