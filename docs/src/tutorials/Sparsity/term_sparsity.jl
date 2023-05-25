# # Term sparsity

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Sparsity/term_sparsity.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Sparsity/term_sparsity.ipynb)
# **Adapted from**: Example 3.5 of [WML20b]
#
# [WML20a] Wang, Jie, Victor Magron, and Jean-Bernard Lasserre.
# *TSSOS: A Moment-SOS hierarchy that exploits term sparsity*.
# arXiv preprint arXiv:1912.08899 (2020).
#
# [WML20b] Wang, Jie, Victor Magron, and Jean-Bernard Lasserre.
# *Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension*.
# arXiv preprint arXiv:2003.03210 (2020).

using Test #src
using DynamicPolynomials
@polyvar x[1:3]

# We would like to find the minimum value of the polynomial

poly = x[1]^2 - 2x[1]*x[2] + 3x[2]^2 - 2x[1]^2*x[2] + 2x[1]^2*x[2]^2 - 2x[2]*x[3] + 6x[3]^2 + 18x[2]^2*x[3] - 54x[2]*x[3]^2 + 142x[2]^2*x[3]^2

# The minimum value of the polynomial can be found to be zero.

import CSDP
solver = CSDP.Optimizer
using SumOfSquares
function sos_min(sparsity)
    model = Model(solver)
    @variable(model, t)
    @objective(model, Max, t)
    con_ref = @constraint(model, poly - t in SOSCone(), sparsity = sparsity)
    optimize!(model)
    return value(t), moment_matrix(con_ref)
end

bound, ν = sos_min(Sparsity.NoPattern())
@test bound ≈ 0 atol=1e-6 #src
bound

# We find the corresponding minimizer `(0, 0, 0)` by matching the moments
# of the moment matrix with a dirac measure centered at this minimizer.

atomic_measure(ν, 1e-6)

# We can see below that the basis contained 6 monomials hence we needed to use 6x6 PSD matrix variables.

@test ν.basis.monomials == [x[1]*x[2], x[2]*x[3], x[1], x[2], x[3], 1] #src
ν.basis

# Using the monomial/term sparsity method of [WML20a] based on cluster completion, we find the same bound.

bound, ν = sos_min(Sparsity.Monomial())
@test bound ≈ 0 atol=1e-6 #src
bound

# Which is not suprising as no sparsity reduction could be performed.

@test length(ν.sub_moment_matrices) == 1 #src
@test ν.sub_moment_matrices[1].basis.monomials == [x[1]*x[2], x[2]*x[3], x[1], x[2], x[3], 1] #src
[sub.basis for sub in ν.sub_moment_matrices]

# Using the monomial/term sparsity method of [WML20b] based on chordal completion, the lower bound is smaller than 0.

bound, ν = sos_min(Sparsity.Monomial(ChordalCompletion()))
@test bound ≈ -0.00355 rtol=1e-3 #src
bound

# However, this bound was obtained with an SDP with 4 matrices of size 3x3.

@test length(ν.sub_moment_matrices) == 4                                  #src
@test ν.sub_moment_matrices[1].basis.monomials == [x[2]*x[3], x[2], x[3]] #src
@test ν.sub_moment_matrices[2].basis.monomials == [x[2]*x[3], x[2], 1]    #src
@test ν.sub_moment_matrices[3].basis.monomials == [x[1], x[2], 1]         #src
@test ν.sub_moment_matrices[4].basis.monomials == [x[1]*x[2], x[1], 1]    #src
[sub.basis for sub in ν.sub_moment_matrices]
