# Sum-of-Squares Programming

This section contains a brief introduction to Sum-of-Squares Programming. For more details, see [BPT12, Las09, Lau09].

## Quadratic forms and Semidefinite programming

The positive semidefiniteness of a ``n \times n`` real symmetric matrix ``Q`` is equivalent to the nonnegativity of the quadratic form ``p(x) = x^\top Q x`` for all vector ``x \in \mathbb{R}^n``.
For instance, the polynomial
```math
x_1^2 + 2x_1x_2 + 5x_2^2 + 4x_2x_3 + x_3^2 = x^\top \begin{pmatrix}1 & 1 & 0\\1 & 5 & 2\\ 0 & 2 & 1\end{pmatrix} x
```
is nonnegative since the matrix of the right-hand side is positive semidefinite.
Moreover, a certificate of nonnegativity can be extracted from the Cholesky decomposition of the matrix:
```math
(x_1 + x_2)^2 + (2x_2 + x_3)^2 = x^\top \begin{pmatrix}1 & 1 & 0\\0 & 2 & 1\end{pmatrix}^\top \begin{pmatrix}1 & 1 & 0\\0 & 2 & 1\end{pmatrix} x
```

## Polynomial nonnegativity and Semidefinite programming

This can be generalized to a polynomial of arbitrary degree.
A polynomial ``p(x)`` is nonnegative is it can be rewritten as ``p(x) = X^\top Q X`` where ``Q`` is a real symmetric positive semidefinite matrix and ``X`` is a vector of monomials.

For instance
```math
x_1^2 + 2x_1^2x_2 + 5x_1^2x_2^2 + 4x_1x_2^2 + x_2^2 = X^\top \begin{pmatrix}1 & 1 & 0\\1 & 5 & 2\\ 0 & 2 & 1\end{pmatrix} X
```
where ``X = (x_1, x_1x_2, x_2)``
Similarly to the previous section, the Cholesky factorization of the matrix can be used to extract a sum of squares certificate of nonnegativity for the polynomial:
```math
(x_1 + x_1x_2)^2 + (2x_1x_2 + x_2)^2 = X^\top \begin{pmatrix}1 & 1 & 0\\0 & 2 & 1\end{pmatrix}^\top \begin{pmatrix}1 & 1 & 0\\0 & 2 & 1\end{pmatrix} X
```

## When is nonnegativity equivalent to sum of squares ?

Determining whether a polynomial is nonnegative is *NP-hard*. The condition of the previous section was only sufficient, that is, there exists nonnegative polynomials that are not sums of squares.
Hilbert showed in 1888 that there are exactly 3 cases for which there is equivalence between the nonnegativity of the polynomials of ``n`` variables and degree ``2d`` and the existence of a sum of squares decomposition.

* ``n = 1`` : Univariate polynomials
* ``2d = 2`` : Quadratic polynomials
* ``n = 2``, ``2d = 4`` : Bivariate quartics

The first explicit example of polynomial that was not a sum of squares was given by Motzkin in 1967:
```math
x_1^4x_2^2 + x_1^2x_2^4 + 1 - 3x_1^2x_2^2 \geq 0 \quad \forall x
```
While it is not a sum of squares, it can still be certified to be nonnegative using sum-of-squares programming by identifying it with a rational sum-of-squares decomposition.
These facts can be verified numerically using this package as detailed in the [motzkin notebook](https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/motzkin.ipynb).

### References

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.

[Las09] Lasserre, J. B.
*Moments, positive polynomials and their applications*
World Scientific, **2009**.

[Lau09] Laurent, M.
*Sums of squares, moment matrices and optimization over polynomials*
Emerging applications of algebraic geometry, Springer, **2009**, 157-270.
