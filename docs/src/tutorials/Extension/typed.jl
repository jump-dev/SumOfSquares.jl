# # Multivariate polynomials implementations

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extension/typed.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extension/typed.ipynb)
# **Contributed by**: Benoît Legat

# The SumOfSquares package is built on top of the [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl)
# abstract interface. [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl/)
# is an implementation of this abstract interface so it can be used with
# SumOfSquares. Moreover, any other implementation can be used as well. To
# illustrate, we solve Examples 3.38 of [BPT12] with
# [TypedPolynomials](https://github.com/JuliaAlgebra/TypedPolynomials.jl),
# another implementation of [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl).
#
# [BPT12] Blekherman, G. & Parrilo, P. A. & Thomas, R. R.
# *Semidefinite Optimization and Convex Algebraic Geometry*.
# Society for Industrial and Applied Mathematics, **2012**.

import TypedPolynomials
TypedPolynomials.@polyvar x y
model = SOSModel(CSDP.Optimizer)
con_ref = @constraint(model, 2x^4 + 5y^4 - x^2*y^2 >= -2(x^3*y + x + 1))
optimize!(model)
solution_summary(model)

# We see that the problem is feasible. The Sum-of-Squares decomposition can be
# obtained as follows:

sos_decomposition(con_ref)

# Why is there several implementations ?
# Depending in the use-case, one implementation may be more appropriate than
# another one. [TypedPolynomials](https://github.com/JuliaAlgebra/TypedPolynomials.jl)
# is faster than [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl/)
# but it requires new compilation whenever the list of variables changes.
# This means that [TypedPolynomials](https://github.com/JuliaAlgebra/TypedPolynomials.jl)
# is not appropriate when the number of variables is dynamic or too large.
# However, for a small number of variables, it can be faster.
# When solving Sum-of-Squares programs, the time is mostly taken by the Semidefinite programming solver.
# The time taken by SumOfSquares/JuMP/MathOptInterface are usually negligible
# or it time is taken by manipulation of JuMP or MathOptInterface functions
# therefore using TypedPolynomials over DynamicPolynomials may not make much difference in most cases.

# One case for which using TypedPolynomials might be adequate is when
# using domain defined by equalities (possibly also with inequalities).
# Indeed, in that case, SumOfSquares computes the corresponding Gröbner basis which
# may take a non-negligible amount of time for large systems of equalities.

# To illustrate this, consider the computation of Gröbner basis for the
# following system from [CLO05, p. 17].
# The time taken by TypedPolynomials is below:
#
# [CLO05] Cox, A. David & Little, John & O'Shea, Donal
# *Using Algebraic Geometry*.
# Graduate Texts in Mathematics, **2005**.
# https://doi.org/10.1007/b138611

using BenchmarkTools
@btime let
    TypedPolynomials.@polyvar x y
    S = @set x^3 * y + x == 2x^2 * y^2 && 3x^4 == y
    SemialgebraicSets.computegröbnerbasis!(S.I)
end

# The time taken by DynamicPolynomials is as follows:

@btime let
    DynamicPolynomials.@polyvar x y
    S = @set x^3 * y + x == 2x^2 * y^2 && 3x^4 == y
    SemialgebraicSets.computegröbnerbasis!(S.I)
end

# We see that TypedPolynomials is faster.
# The time is still negligible for this small system but for larger systems, choosing TypedPolynomials may be helpful.
# We can use this system in a Sum-of-Squares constraint as follows:

TypedPolynomials.@polyvar x y
S = @set x^3 * y + x == 2x^2 * y^2 && 3x^4 == y
poly = -6x - 4y^3 + 2x*y^2 + 6x^3 - 3y^4 + 13x^2 * y^2
model = Model(CSDP.Optimizer)
con_ref = @constraint(model, poly in SOSCone(), domain = S)
optimize!(model)
solution_summary(model)

# We obtain the following decomposition:

dec = sos_decomposition(con_ref, 1e-3)

# We can verify that it is correct as follows:

# TODO remove `polynomial` #src
@test isapproxzero(rem(polynomial(dec) - poly, S.I), ztol = 1e-5) #src
rem(polynomial(dec) - poly, S.I)
