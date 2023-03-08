using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra
using SCS
using JuMP

function test(solver, d, obj, I, lagrangian_monos)
    model = SOSModel(solver)
    @variable(model, t)
    if lagrangian_monos isa Float64
		display(I)
        c = @constraint(model, t >= obj, domain = I, maxdegree = d)
		@show moi_set(constraint_object(c)).certificate
		@show moi_set(constraint_object(c)).domain
		@show I.I.gröbnerbasis
		display(I)
    else
        p = equalities(I)
        @variable(model, multipliers[eachindex(p)], Poly(lagrangian_monos))
        @constraint(model, t >= obj + multipliers ⋅ p)
    end
    @objective(model, Min, t)
    optimize!(model)
    solution_summary(model)
end

function sol(solver, d, ztol)
    X = [0.5459627556242905, 1.7288950489507429, 0.7167681447476535]

    Y = [-54.06002080721971, 173.77393714162503, 71.48154370522498]

    n = 3

    @polyvar W[1:n] α β

	I = @set(sum([W[i] * X[i] * (β * X[i] + α - Y[i]) for i = 1:n]) == 0, library = Buchberger(ztol))
    #I = @set sum([W[i] * X[i] * (β * X[i] + α - Y[i]) for i = 1:n]) == 0
	for i in 1:n
		addequality!(I, W[i] - W[i]^2)
	end

    obj = sum(W)

	if ztol === nothing
		monos = monomials([α, β], 0:d)
	else
		monos = ztol
	end

	test(solver, d, obj, I, monos)
end

function test(solver)
    X = [0.5459627556242905, 1.7288950489507429, 0.7167681447476535]

    Y = [-54.06002080721971, 173.77393714162503, 71.48154370522498]

    model = SOSModel(solver)

    n = 3

    @polyvar W[1:n] α β

    d = 2

    ab_mons = monomials([α, β], 0:d)

    # build boolean ideal
    @variable(model, p_coeff[1:n], Poly(ab_mons))
    p = sum([p_coeff[i] * (W[i] - W[i]^2) for i = 1:n])

    # build objective function
    @variable(model, t)
    obj = sum(W)

    # one more ideal generator
    r1_gen = sum([W[i] * X[i] * (β * X[i] + α - Y[i]) for i = 1:n])

    # multiplier for the ideal
    @variable(model, r1_coeff, Poly(ab_mons))
    r1 = r1_coeff * r1_gen

    @constraint(model, t - obj - p - r1 >= 0)
    @objective(model, Min, t)

    optimize!(model)

    solution_summary(model)
end
