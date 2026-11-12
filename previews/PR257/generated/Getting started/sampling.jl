using DynamicPolynomials
using SumOfSquares
import MultivariateBases as MB

using Dualization
import Hypatia
import SCS

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 3

model = Model(dual_optimizer(SCS.Optimizer))
set_silent(model)
@variable(model, γ)
@objective(model, Max, γ)
@constraint(model, p - γ in SOSCone(), zero_basis = BoxSampling([-1], [1]))
optimize!(model)
solution_summary(model)

print_active_bridges(model)

set_optimizer(model, dual_optimizer(Hypatia.Optimizer))
optimize!(model)
solution_summary(model)

print_active_bridges(model)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
