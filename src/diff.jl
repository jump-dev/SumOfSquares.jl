# I do not use it but I import the function to add a method
import Calculus.differentiate

function differentiate(p::VecPolynomial, x::PolyVar)
  # grlex order preserved
  i = findfirst(vars(p), x)
  if i == 0
    0
  else
    keep = find([z[i] > 0 for z in p.x.Z])
    Z = [copy(p.x.Z[i]) for i in keep]
    a = p.a[keep] # does a copy
    for j in 1:length(Z)
      a[j] *= Z[j][i]
      Z[j][i] -= 1
    end
    VecPolynomial(a, MonomialVector(vars(p), Z))
  end
end

function differentiate(p::VecPolynomial, xs::Vector{PolyVar})
  [differentiate(p, x) for x in xs]
end

differentiate(p::MatPolynomial, x) = differentiate(VecPolynomial(p), x)
