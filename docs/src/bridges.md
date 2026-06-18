# Bridging mechanism

[Bridges](https://jump.dev/MathOptInterface.jl/stable/submodules/Bridges/overview/)
are MathOptInterface's mechanism for taking a constraint or constrained
variable expressed in one set and rewriting it in terms of other sets, until
something that the underlying solver supports natively is reached.

In SumOfSquares, the user-facing object — a JuMP `@constraint` such as
```julia
@constraint(model, p in SOSCone(), domain = (@set x^2 + y^2 == 1 && x >= 0))
```
is several bridge hops away from the cone that the solver will actually see.
This page walks through those hops, so that you can predict — and tune — the
chain of reformulations between the SOS constraint that you wrote and the
PSD (or rotated-second-order-cone, or rank-1) constraint that hits the
solver.

The high-level shape of the chain is:

```text
SOSPolynomialSet ────► WeightedSOSCone ────► (PSD-like cone) ────► solver
                ▲                       ▲
            constraint               variable
              bridges                  bridges
```

The first hop is always a *constraint bridge*. The second hop, which produces
the gram-matrix PSD constraints, can be either a *variable bridge* (e.g.
[`SumOfSquares.Bridges.Variable.KernelBridge`](@ref)) or a *constraint
bridge* (e.g. [`SumOfSquares.Bridges.Constraint.ImageBridge`](@ref)). The
choice is driven by which form the underlying solver supports more directly,
and the bridge graph picks the cheapest path.

## From a JuMP constraint to a `WeightedSOSCone`

[`SumOfSquares.SOSPolynomialSet`](@ref) is the MOI set that JuMP creates
when you ask for an SOS constraint over a domain. The domain is carried as
a parameter:

  * `SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}` — equality
    constraints only (or no `domain` at all);
  * `SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet}` — at least
    one inequality constraint.

Two bridges turn this into a [`SumOfSquares.WeightedSOSCone`](@ref):

  * [`SumOfSquares.Bridges.Constraint.SOSPolynomialBridge`](@ref) handles the
    algebraic case. It reduces the polynomial modulo the algebraic ideal,
    asks the certificate for a single gram basis, and produces a
    `WeightedSOSCone` with weight ``1``.
  * [`SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge`](@ref)
    handles the basic-semialgebraic case (Putinar-style). It collects one
    multiplier basis per inequality ``g_i`` and produces a single
    `WeightedSOSCone` with weights ``[1, g_1, g_2, \ldots]`` and gram bases
    ``[\sigma_0, \sigma_1, \sigma_2, \ldots]``. The σ_0 basis comes from the
    inner ideal certificate.

After this hop, every SOS constraint is a single `WeightedSOSCone`
constraint, regardless of whether the user supplied a domain.

## From a `WeightedSOSCone` to PSD

There are two natural ways to lower a `WeightedSOSCone`:

  * **Variable side** — add the gram matrices ``Q_i`` as PSD-constrained
    variables and express each polynomial coefficient as an affine function
    of the ``Q_i``. This is what
    [`SumOfSquares.Bridges.Variable.KernelBridge`](@ref) does (and what
    [`SumOfSquares.Bridges.Variable.LowRankBridge`](@ref) does for the
    rank-1 case, see below).

  * **Constraint side** — pick anchor entries of each ``Q_i`` so that they
    are *fixed by* the polynomial coefficients (no new variables), introduce
    slack variables only when two gram entries produce the same monomial,
    and finally add the matrices as PSD constraints. This is what
    [`SumOfSquares.Bridges.Constraint.ImageBridge`](@ref) does; the
    docstring contains a worked example.

Mathematically, if ``Σ`` is the SOS cone and ``S`` the PSD cone, then
``Σ = A(S)`` for a linear map ``A``. The variable bridge gives the user
``A(q)`` directly; the constraint bridge gives the user the preimage
``q \in A^{-1}(p) = A^\dagger p + \ker(A)`` instead.

Which one is cheaper depends on the solver:

| Solver style                                                                                 | Wins |
| -------------------------------------------------------------------------------------------- | ---- |
| Native PSD-constrained variables (e.g. Mosek's `barvar`, Hypatia)                            | `KernelBridge` |
| PSD as a constraint only (e.g. Clarabel, SCS)                                                | `ImageBridge` |
| Specialised low-rank PSD format (e.g. Loraine, Burer–Monteiro)                               | `LowRankBridge` + LowRankOpt bridges |

The bridge graph uses `MOI.Bridges.bridging_cost` to pick. For the SOS
pipeline, `ImageBridge` carries an explicit `bridging_cost(::Type{<:ImageBridge}) = 2.0`
so that `KernelBridge` wins for solvers that support constrained PSD
variables (and vice versa for solvers that only support PSD as a
constraint).

## How Hypatia gets specialised cones

Hypatia natively supports `MOI.PositiveSemidefiniteConeTriangle` *and* a
range of more specialised cones (`HypatiaSOSCone`, etc.). Because the SOS
pipeline lands in `WeightedSOSCone` after the first hop,
[`SumOfSquares.Bridges.Variable.KernelBridge`](@ref) wins (Hypatia accepts
PSD as constrained variables) and the chain stops there:

```text
SOSPolynomialSet
   └─SOSPolynomialBridge──► WeightedSOSCone
                                 └─KernelBridge──► PositiveSemidefiniteConeTriangle (Hypatia native)
```

When the cone is structurally special (Lagrange basis with rank-1 factors,
see next section), the chain takes a different exit through
`LowRankBridge` and Hypatia still ends up receiving a cone it understands
natively — the bridge graph is what makes this transparent.

## How LowRankOpt exploits rank-1 structure

When the gram basis is a [`MultivariateBases.LagrangeBasis`](@ref) — the
polynomial is sampled at a fixed set of points — the relation between
gram-matrix entries and polynomial coefficients factorises through a list
of rank-1 outer products, one per Lagrange node:

```math
p_j = \sum_i w_i(x_j) \langle u_{i,j} u_{i,j}^\top, Q_i \rangle
```

where ``u_{i,j}`` is a column of the basis transformation. The
[LowRankOpt](https://github.com/blegat/LowRankOpt.jl) package provides MOI
sets that carry this rank-1 structure all the way down to the solver:

  * `LRO.SetDotProducts{W,S,V}` — the set of vectors
    ``(\langle a_1, x \rangle, \ldots, \langle a_m, x \rangle)`` for ``x``
    in some PSD-like set `S`. The factors `a_k` are exposed as
    `LRO.TriangleVectorization{LRO.Factorization{...}}`, i.e. rank-1
    matrices given by their factor and weight rather than by their full
    matrix.
  * `LRO.LinearCombinationInSet{W,S,V}` — the dual side: vectors ``y``
    such that ``\sum_k y_k a_k`` lies in ``S``.

[`SumOfSquares.Bridges.Variable.LowRankBridge`](@ref) reformulates a
`WeightedSOSCone` whose gram basis is a `LagrangeBasis` directly into
`LRO.SetDotProducts{WITHOUT_SET}` wrapping the appropriate
`SOS.matrix_cone`, with each Lagrange node giving rise to one rank-1
factor. So the chain looks like:

```text
SOSPolynomialSet
   └─SOSPolynomial(InSemialgebraicSet)Bridge──► WeightedSOSCone (Lagrange)
                                                       └─LowRankBridge──► SetDotProducts{WITHOUT_SET, PSD, Factorization}
```

From there, what happens next is solver-dependent.

### Solvers that natively support the rank-1 format

* **Loraine.jl** — implements an interior-point method whose Hessian
  computations exploit the rank-1 structure of the constraint matrices.
  It supports `LRO.SetDotProducts` directly, so the chain stops at the
  output of `LowRankBridge`.

* **LowRankOpt.BurerMonteiro** — a Burer–Monteiro-style nonlinear-program
  back-end. It also consumes `LRO.SetDotProducts` directly.

### Solvers that don't

If the solver only knows ordinary PSD, the LowRankOpt bridge layer takes
over and unfolds the structure step by step. The most relevant LowRankOpt
bridges, in roughly the order they fire, are:

  * `LowRankOpt.Bridges.Variable.AppendSetBridge` — converts between
    `SetDotProducts{WITHOUT_SET}` and `SetDotProducts{WITH_SET}` (the
    `WITH_SET` form additionally carries the underlying PSD slice as
    explicit variables).
  * `LowRankOpt.Bridges.Variable.DotProductsBridge` — drops the rank-1
    factorisation and rewrites the dot products as explicit affine
    combinations of PSD variables. This is the bridge that "loses" the
    rank-1 information for solvers that can't use it.
  * `LowRankOpt.Bridges.Variable.ToPositiveBridge` and
    `LowRankOpt.Bridges.Variable.ToRankOneBridge` — convert between
    different concrete `LRO.Factorization` types (e.g. rank-`r` factors
    expressed as `r` rank-1 factors).
  * `LowRankOpt.Bridges.Constraint.LinearCombinationBridge` and
    `LowRankOpt.Bridges.Constraint.AppendZeroBridge` — the constraint-side
    duals of the variable bridges above, used when the underlying solver
    accepts the cones only as constraints.

The net effect is that the SOS chain stays the same up to the
`LRO.SetDotProducts` level, and a fallback to plain
`MOI.PositiveSemidefiniteConeTriangle` only happens at the very end if no
solver in the stack can do better.

## Choosing a different bridge

The bridge graph is driven entirely by the *types* the solver supports and
by `MOI.Bridges.bridging_cost`. There are two practical ways to nudge the
chain:

  * Use `JuMP.remove_bridge(model, SomeBridge)` to disable a bridge — for
    example, `JuMP.remove_bridge(model, SumOfSquares.Bridges.Variable.KernelBridge{Float64})`
    forces the chain to go through `ImageBridge` even for solvers that
    accept PSD-constrained variables. This is occasionally useful for
    debugging or for benchmarking.
  * Define a new bridge with a lower cost than the existing one. This is
    how solver packages opt into specialised cones: by providing a
    `MOI.supports_constrained_variable` (or `MOI.supports_constraint`) for
    the special set and a bridge with a cost that makes the bridge graph
    prefer the specialised route over the generic one.

In both cases, you can inspect the resulting chain with
`MOI.Bridges.print_active_bridges(JuMP.backend(model))` to confirm the
exact path that a given constraint takes.
