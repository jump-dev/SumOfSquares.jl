import Base.dot, Base.(.^)

(.^)(p::TermType, i::Int) = p^i

function (*)(x::PolyVar, y::PolyVar)
  if x === y
    Monomial([x], [2])
  else
    Monomial([x,y], [1,1])
  end
end
function (*)(x::PolyVar, y::Monomial)
  i = 1
  vars = y.vars
  n = length(vars)
  while i <= n && x > vars[i]
    i += 1
  end
  if i > n
    Monomial([vars; x], [y.z; 1])
  elseif x == vars[i]
    znew = copy(y.z)
    znew[i] += 1
    Monomial(vars, znew)
  else
    znew = copy(y.z)
    znew[i] += 1
    Monomial([vars[1:i-1]; x; vars[i:end]], [y.z[1:i-1]; 1; y.z[i:end]])
  end
end
(*)(x::Monomial, y::PolyVar) = y * x
function (*)(x::Monomial, y::Monomial)
  if x.vars == y.vars
    Monomial(x.vars, x.z + y.z)
  else
    if x.vars == y.vars
      Monomial(x.vars, x.z+y.z)
    else
      allvars, maps = myunion([x.vars, y.vars])
      z = zeros(Int, length(allvars))
      z[maps[1]] += x.z
      z[maps[2]] += y.z
      Monomial(allvars, z)
    end
  end
end

*(α, x::PolyVar) = Term(α, Monomial(x))
*(α, x::Monomial) = Term(α, x)
*(α, x::Term) = Term(T(α)*x.α, x.x)
*{T<:Union{PolyVar,Monomial,Term}}(x::T, α) = α * x
*(x::Term, y::Term) = Term(x.α*y.α, x.x*y.x)
*(x::Term, y::PolyVar) = Term(x.α, x.x*y)
*(x::PolyVar, y::Term) = y * x
*(x::Term, y::Monomial) = Term(x.α, x.x*y)
*(x::Monomial, y::Term) = y * x

*(α, p::MatPolynomial) = α * VecPolynomial(p)
*(p::MatPolynomial, α) = VecPolynomial(p) * α

*(α, p::VecPolynomial) = VecPolynomial(α*p.a, p.x)

function *(p::VecPolynomial, q::VecPolynomial)
  if iszero(p)
    zero(q)
  elseif iszero(q)
    zero(p)
  else
    samevars = vars(p) == vars(q)
    if samevars
      allvars = vars(p)
    else
      allvars, maps = myunion([vars(p), vars(q)])
    end
    N = length(p)*length(q)
    Z = Vector{Vector{Int}}(N)
    T = typeof(p.a[1]*q.a[1])
    a = Vector{T}(N)
    i = 0
    for u in p
      for v in q
        if samevars
          z = u.x.z + v.x.z
        else
          z = zeros(Int, length(allvars))
          z[maps[1]] += u.x.z
          z[maps[2]] += v.x.z
        end
        i += 1
        Z[i] = z
        a[i] = u.α * v.α
      end
    end
    vecpolynomialclean(allvars, a, Z)
  end
end

dot(p::PolyType, q::PolyType) = p * q
dot(α, p::PolyType) = α * p
dot(p::PolyType, α) = p * α

myminivect{T}(x::T, y::T) = [x, y]
function myminivect{S,T}(x::S, y::T)
  U = promote_type(S, T)
  [U(x), U(y)]
end

function (+)(x::Term, y::Term)
  if x.x == y.x
    Term(x.α+y.α, x.x)
  else
    VecPolynomial(myminivect(x.α,y.α), [x.x,y.x])
  end
end

function (-)(x::Term, y::Term)
  if x.x == y.x
    Term(x.α-y.α, x.x)
  else
    VecPolynomial(myminivect(x.α,-y.α), [x.x,y.x])
  end
end

(+){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) + Term(y)
(-){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) - Term(y)

function plusorminus{S,T}(x::TermContainer{S}, y::TermContainer{T}, isplus)
  varsvec = [vars(x), vars(y)]
  allvars, maps = myunion(varsvec)
  nvars = length(allvars)
  U = promote_type(S, T)
  a = Vector{U}()
  Z = Vector{Vector{Int}}()
  for (i, tc) in enumerate([x,y])
    for (j, t) in enumerate(tc)
      added = false
      z = zeros(Int, nvars)
      z[maps[i]] = t.x.z
      if i == 1 || isplus
        α = t.α
      else
        α = -t.α
      end
      for k in 1:length(Z)
        if Z[k] == z
          a[k] += α
          added = true
          break
        end
      end
      if !added
        push!(a, α)
        push!(Z, z)
      end
    end
  end
  VecPolynomial(a, MonomialVector(allvars, Z))
end


(+)(x::TermContainer, y::TermContainer) = plusorminus(x, y, true)
(-)(x::TermContainer, y::TermContainer) = plusorminus(x, y, false)
(+){S<:Union{Monomial,PolyVar},T}(x::TermContainer{T}, y::S) = x + Term{T}(y)
(+){S<:Union{Monomial,PolyVar},T}(x::S, y::TermContainer{T}) = Term{T}(x) + y

(+)(x::TermContainer, y::MatPolynomial) = x + VecPolynomial(y)
(+)(x::MatPolynomial, y::TermContainer) = VecPolynomial(x) + y
(-)(x::TermContainer, y::MatPolynomial) = x - VecPolynomial(y)
(-)(x::MatPolynomial, y::TermContainer) = VecPolynomial(x) - y

iszero{T}(x::T) = x == zero(T)
iszero(t::Term) = iszero(t.α)
iszero(p::VecPolynomial) = isempty(p.x)
iszero(p::MatPolynomial) = isempty(p.x)

(-){S<:Union{Monomial,PolyVar},T}(x::TermContainer{T}, y::S) = x - Term{T}(y)
(-){S<:Union{Monomial,PolyVar},T}(x::S, y::TermContainer{T}) = Term{T}(x) - y

# Avoid adding a zero constant that might artificially increase the Newton polytope
(+)(x::PolyType, y) = iszero(y) ? x : x + Term(y)
(+)(x, y::PolyType) = iszero(x) ? y : Term(x) + y
(-)(x::PolyType, y) = iszero(y) ? x : x - Term(y)
(-)(x, y::PolyType) = iszero(x) ? -y : Term(x) - y

(-)(x::PolyVar) = Term(-1, Monomial(x))
(-)(x::Monomial) = Term(-1, x)
