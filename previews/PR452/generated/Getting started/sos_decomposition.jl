using DynamicPolynomials
using SumOfSquares
import CSDP

@polyvar x y
p = x^2 - x*y^2 + y^4 + 1

model = SOSModel(CSDP.Optimizer)
@constraint(model, cref, p >= 0)

optimize!(model)

sos_dec = sos_decomposition(cref, 1e-4)

polynomial(sos_dec, Float32)

gram = gram_matrix(cref)

keys_as_monomials(gram.basis)' * gram.Q * keys_as_monomials(gram.basis)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
