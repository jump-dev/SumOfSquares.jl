using LinearAlgebra
using DynamicPolynomials
@polyvar s[1:3]

quad = 2 * sum(s .^ 2) + 8 * (s[1] * s[2] + s[1] * s[3] + s[2] * s[3])

using SumOfSquares
import Clarabel
solver = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

function sos_lower_bound(p; multiplier = 1)
    model = SOSModel(solver)
    @variable(model, γ)
    @objective(model, Max, γ)
    con_ref = @constraint(model, multiplier * (p - γ) >= 0)
    optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL
    return value(γ), con_ref
end

m_A = quad + sum(s .^ 4)

γ_A, _ = sos_lower_bound(m_A)
γ_A

γ_A_mult, con_ref = sos_lower_bound(m_A, multiplier = sum(s .^ 2))
γ_A_mult

atoms = atomic_measure(moment_matrix(con_ref), 1e-4)

minimizers = [atom.center for atom in atoms.atoms if atom.weight > 1e-4]

m_A_min = minimum(x -> m_A(s => x), minimizers)
m_A_min

for σ in [1, 10, 100]
    m = quad + σ / 4 * sum(s .^ 4)
    γ, _ = sos_lower_bound(m)
    println("σ = $σ: lower bound = $γ, minimum = $(4m_A_min / σ)")
end

σ = 4
m_E = quad + σ / 4 * sum(s .^ 2)^2

γ_E, con_ref = sos_lower_bound(m_E)
γ_E

atoms = atomic_measure(moment_matrix(con_ref), 1e-4)

for atom in atoms.atoms
    x = atom.center
    println("sum = $(sum(x)), norm = $(norm(x))")
end

cubic = 0.1 * sum(s .^ 3)
m_A2 = m_A + cubic
γ_A2, _ = sos_lower_bound(m_A2)
γ_A2

γ_A2_mult, con_ref = sos_lower_bound(m_A2, multiplier = sum(s .^ 2))
γ_A2_mult

for σ in [0.1, 1, 4, 10, 40]
    m = quad + cubic + σ / 4 * sum(s .^ 2)^2
    γ, c = sos_lower_bound(m)
    μ = atomic_measure(moment_matrix(c), 1e-4)
    m_min = minimum(atom -> m(s => atom.center), μ.atoms)
    println("σ = $σ: lower bound = $γ, value at extracted minimizers = $m_min")
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
