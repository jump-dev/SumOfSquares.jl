export Term, VecPolynomial, MatPolynomial

abstract TermContainer{T}

type Term{T} <: TermContainer{T}
  α::T
  x::Monomial
end
Term(t::Term) = t
Term(x::Monomial) = Term{Int}(x)
Term(x::PolyVar) = Term(Monomial(x))
Term{T}(α::T) = Term{T}(α, Monomial())
(::Type{Term{T}}){T}(t::Term{T}) = t
(::Type{Term{T}}){T}(t::Term) = Term{T}(T(t.α), t.x)
(::Type{Term{T}}){T}(x::Monomial) = Term{T}(one(T), x)
(::Type{Term{T}}){T}(x::PolyVar) = Term{T}(Monomial(x))
(::Type{Term{T}}){T}(α) = Term(T(α))

vars(t::Term) = vars(t.x)

length(::Term) = 1
isempty(::Term) = false
start(::Term) = false
done(::Term, state) = state
next(x::Term, state) = (x, true)

type VecPolynomial{T} <: TermContainer{T}
  a::Vector{T}
  x::MonomialVector

  function VecPolynomial(a::Vector{T}, x::MonomialVector)
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
VecPolynomial{T}(a::Vector{T}, x::MonomialVector) = VecPolynomial{T}(a, x)

function vecpolynomialclean{T}(vars::Vector{PolyVar}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
  σ = sortperm(Zdup)
  Z = Vector{Vector{Int}}()
  a = Vector{T}()
  i = 0
  j = 1
  while j <= length(adup)
    k = σ[j]
    if j == 1 || Zdup[k] != Zdup[σ[j-1]]
      push!(Z, Zdup[k])
      push!(a, adup[k])
      i += 1
    else
      a[i] += adup[k]
    end
    j += 1
  end
  VecPolynomial(a, MonomialVector(vars, Z))
end

vars(p::VecPolynomial) = vars(p.x)

length(x::VecPolynomial) = length(x.a)
isempty(x::VecPolynomial) = length(x) > 0
start(::VecPolynomial) = 1
done(x::VecPolynomial, state) = length(x) < state
next(x::VecPolynomial, state) = (Term(x.a[state], x.x[state]), state+1)

type MatPolynomial{T}
  Q::Vector{T}
  x::MonomialVector
end

function trimap(i, j, n)
  div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end

function (::Type{MatPolynomial{T}}){T}(f::Function, x::MonomialVector)
  n = length(x)
  Q = Vector{T}(trimap(n, n, n))
  for i in 1:n
    for j in i:n
      Q[trimap(i,j,n)] = f(i,j)
    end
  end
  MatPolynomial{T}(Q, x)
end

function getindex(p::MatPolynomial, I::NTuple{2,Int})
  i, j = I
  if i < j
    i, j = (j, i)
  end
  n = length(p.x)
  p.Q[trimap(i,j,n)]
end

function VecPolynomial{T}(p::MatPolynomial{T})
  n = length(p.x)
  N = trimap(n, n, n)
  Z = Vector{Vector{Int}}(N)
  U = typeof(2*p.Q[1] + p.Q[1])
  a = Vector{U}(N)
  for i in 1:n
    for j in i:n
      k = trimap(i, j, n)
      Z[k] = p.x.Z[i] + p.x.Z[j]
      if i == j
        a[k] = p.Q[k]
      else
        a[k] = 2*p.Q[k]
      end
    end
  end
  vecpolynomialclean(p.x.vars, a, Z)
end
