using DynamicPolynomials
using SumOfSquares
import CSDP

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 3

model = SOSModel(CSDP.Optimizer)
@variable(model, σ)
@constraint(model, cref, p >= σ)
@objective(model, Max, σ)

optimize!(model)
solution_summary(model)

sos_dec = sos_decomposition(cref, 1e-4)

ν = moment_matrix(cref)
ν.Q

η = extractatoms(ν, 1e-4)
minimizers = [η.atoms[1].center; η.atoms[2].center]

η1 = moment_matrix(dirac(monomials(x, 0:4), x => round(minimizers[1])), ν.basis.monomials)
η1.Q

η2 = moment_matrix(dirac(monomials(x, 0:4), x => round(minimizers[2])), ν.basis.monomials)
η2.Q

Q12 = η1.Q * η.atoms[1].weight + η2.Q * η.atoms[2].weight

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

