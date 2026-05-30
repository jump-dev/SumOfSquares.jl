# BMKKTSystem — matrix-free MadNLP KKT system for LowRankOpt.BurerMonteiro.Model
# Routes through `BurerMonteiro.Model`'s `NLPModels.hprod!`/`jprod!`/`jtprod!`
# (so the Stage-3 batched-FFT path stays live) and solves each Newton step
# with `Krylov.minres!` (no preconditioner).
#
# Closely modeled after `CompressedSensingIPM/src/fft_kkt.jl`.

import MadNLP
import Krylov
import NLPModels
import LinearAlgebra
import LowRankOpt as LRO

mutable struct BMKKTSystem{QLP,T,VT,NLP,LS} <:
               MadNLP.AbstractReducedKKTSystem{T,VT,Matrix{T},MadNLP.ExactHessian{T,VT}}
    nlp::NLP
    n::Int
    m::Int
    # MadNLP-standard diagonals it updates between iterations
    reg::VT
    pr_diag::VT
    du_diag::VT
    l_diag::VT
    u_diag::VT
    l_lower::VT
    u_lower::VT
    ind_lb::Vector{Int}
    ind_ub::Vector{Int}
    # Current Lagrangian multipliers — stashed at `eval_lag_hess_wrapper!` time
    # and read by `mul!` (the BM Hessian depends on `y` but not on `x`).
    current_y::VT
    # Current primal iterate — stashed at `eval_jac_wrapper!` time. The BM
    # constraint c(x) is *quadratic* in `x` (since x is the rank factor /
    # square-scalar pre-image), so its Jacobian depends on the linearization
    # point and we must pass the actual IPM iterate to `jprod!`/`jtprod!`,
    # not the Krylov direction.
    current_x::VT
    # Krylov state
    linear_solver::LS
    krylov_iterations::Vector{Int}
    krylov_residuals::Vector{Float64}
    # Buffers reused inside `mul!`
    hv_buf::VT
    jv_buf::VT
    jtv_buf::VT
    # Diagonal preconditioner for Krylov. Length `n+m`. Refreshed before
    # each `krylov_solve!` from `reg + pr_diag + |diag(H)|` on the primal
    # block and `max(|du_diag|, 1)` on the dual block. `diag(H)` is
    # extracted by probing `hprod!(e_i)`; cheap for our smoke-test sizes
    # but `O(n)` FFTs per IPM iter for the larger benchmark.
    precond_diag::VT
    precond_probe::VT      # buffer for the probe vector
end
# Wrapper so Krylov treats `precond_diag` as a left preconditioner via
# `LinearAlgebra.ldiv!` (we set `ldiv = true` in `solve_kkt!`).
struct BMJacobiPrecond{VT} <: AbstractMatrix{Float64}
    diag::VT
end
Base.size(M::BMJacobiPrecond) = (length(M.diag), length(M.diag))
Base.eltype(::BMJacobiPrecond{VT}) where {VT} = eltype(VT)
LinearAlgebra.ldiv!(y::AbstractVector, M::BMJacobiPrecond, x::AbstractVector) =
    (y .= x ./ M.diag; return y)
LinearAlgebra.ldiv!(M::BMJacobiPrecond, x::AbstractVector) = (x ./= M.diag; return x)

# `qlp=true` → `Krylov.MinresQlpWorkspace` (robust on singular K, default).
# `qlp=false` → `Krylov.MinresWorkspace` (cheaper per iter, but falls back to
# the min-norm least-squares solution when K is rank-deficient — yields
# bogus Newton directions there).
function BMKKTSystem(nlp; qlp::Bool = true, T = Float64, VT = Vector{T})
    n = NLPModels.get_nvar(nlp)
    m = NLPModels.get_ncon(nlp)
    workspace = qlp ?
        Krylov.MinresQlpWorkspace(n + m, n + m, VT) :
        Krylov.MinresWorkspace(n + m, n + m, VT)
    return BMKKTSystem{qlp,T,VT,typeof(nlp),typeof(workspace)}(
        nlp, n, m,
        VT(undef, n),                              # reg
        VT(undef, n),                              # pr_diag
        VT(undef, m),                              # du_diag
        VT(undef, 0),                              # l_diag (no lb)
        VT(undef, 0),                              # u_diag (no ub)
        VT(undef, 0),                              # l_lower
        VT(undef, 0),                              # u_lower
        Int[],                                     # ind_lb
        Int[],                                     # ind_ub
        zeros(T, m),                               # current_y
        zeros(T, n),                               # current_x
        workspace,
        Int[],
        Float64[],
        VT(undef, n),                              # hv_buf
        VT(undef, m),                              # jv_buf
        VT(undef, n),                              # jtv_buf
        ones(T, n + m),                            # precond_diag
        zeros(T, n),                               # precond_probe
    )
end

function MadNLP.create_kkt_system(
    ::Type{BMKKTSystem{QLP}},
    cb::MadNLP.AbstractCallback{T,VT},
    linear_solver::Type;
    opt_linear_solver = MadNLP.default_options(linear_solver),
    hessian_approximation = MadNLP.ExactHessian,
    qn_options = MadNLP.QuasiNewtonOptions(),
) where {QLP,T,VT}
    # Only supported with `square_scalars=true` (no bounds, equality
    # constraints only). Verify here so misconfiguration surfaces early.
    @assert isempty(cb.ind_ineq) "BMKKTSystem assumes equality-only constraints"
    @assert isempty(cb.ind_lb)   "BMKKTSystem assumes no lower-bound variables (square_scalars=true)"
    @assert isempty(cb.ind_ub)   "BMKKTSystem assumes no upper-bound variables (square_scalars=true)"
    return BMKKTSystem(cb.nlp; qlp = QLP, T = T, VT = VT)
end
# Convenience dispatch — `kkt_system = BMKKTSystem` (no type param) ↔ QLP=true.
MadNLP.create_kkt_system(::Type{BMKKTSystem}, cb, ls; kws...) =
    MadNLP.create_kkt_system(BMKKTSystem{true}, cb, ls; kws...)

MadNLP.num_variables(kkt::BMKKTSystem) = kkt.n
MadNLP.get_hessian(::BMKKTSystem) = nothing
MadNLP.get_jacobian(::BMKKTSystem) = nothing

# Pretend the Krylov workspace is a "factorization" with no inertia info.
const _KrylovWS = Union{Krylov.MinresWorkspace,Krylov.MinresQlpWorkspace}
MadNLP.is_inertia(::_KrylovWS) = true
MadNLP.inertia(::_KrylovWS) = (0, 0, 0)
MadNLP.introduce(::Krylov.MinresQlpWorkspace) = "Krylov.MINRES-QLP"
MadNLP.introduce(::Krylov.MinresWorkspace)    = "Krylov.MINRES"
MadNLP.improve!(::_KrylovWS) = true
MadNLP.factorize!(::_KrylovWS) = nothing
MadNLP.is_inertia_correct(::BMKKTSystem, _, _, _) = true

Base.eltype(::BMKKTSystem{QLP,T}) where {QLP,T} = T
Base.size(kkt::BMKKTSystem) = (kkt.n + kkt.m, kkt.n + kkt.m)
Base.size(kkt::BMKKTSystem, ::Int) = kkt.n + kkt.m

function MadNLP.initialize!(kkt::BMKKTSystem{QLP,T}) where {QLP,T}
    fill!(kkt.reg, one(T))
    fill!(kkt.pr_diag, one(T))
    fill!(kkt.du_diag, zero(T))
    fill!(kkt.current_y, zero(T))
    fill!(kkt.current_x, zero(T))
    return
end

# Never assemble Jacobian or Hessian — but DO stash the current iterate so
# our matrix-free `mul!` linearizes the (nonlinear) BM constraint at it.
function MadNLP.eval_jac_wrapper!(
    ::MadNLP.MadNLPSolver, kkt::BMKKTSystem, x::MadNLP.PrimalVector,
)
    copyto!(kkt.current_x, MadNLP.full(x))
    return
end
function MadNLP.eval_lag_hess_wrapper!(
    ::MadNLP.MadNLPSolver,
    kkt::BMKKTSystem,
    ::MadNLP.PrimalVector,
    l::AbstractVector;
    is_resto = false,
)
    # MadNLP convention: `L(x, y) = obj_weight·f(x) + y'·c(x)` (dual feasibility
    # `∇f + Jᵀy = 0`, see `IPM/kernels.jl:247`).
    # NLPModels convention: `L(x, y) = obj_weight·f(x) − y'·c(x)`, so its
    # `hprod!` returns `(obj_weight·∇²f − Σ y_i ∇²c_i)·v`.
    # To get MadNLP's Hessian `∇²f + Σ y_i ∇²c_i` we must pass `−l`.
    kkt.current_y .= .-l
    return
end

MadNLP.compress_jacobian!(::BMKKTSystem) = nothing
MadNLP.compress_hessian!(::BMKKTSystem) = nothing
MadNLP.build_kkt!(::BMKKTSystem) = nothing
function MadNLP.factorize_wrapper!(
    solver::MadNLP.MadNLPSolver{T,VT,IT,KKT},
) where {T,VT,IT,KKT<:BMKKTSystem}
    MadNLP.build_kkt!(solver.kkt)
    return true
end

# Augmented KKT matvec on the `[Δx; Δy]` layout:
#   yp = β·yp + α·(H·xp + Aᵀ·xd + pr_diag·xp)
#   yd = β·yd + α·(A·xp − du_diag·xd)
function _kkt_apply!(
    yp::AbstractVector, yd::AbstractVector,
    kkt::BMKKTSystem,
    xp::AbstractVector, xd::AbstractVector,
    alpha::Number, beta::Number,
)
    # Augmented KKT operator (no bounds, no slacks; reg + du_diag added by
    # the trailing `_kktmul!`, which mirrors MadNLP's standard sparse path):
    #   yp = β·yp + α·(H·xp + Aᵀ·xd)
    #   yd = β·yd + α·(A·xp)
    # Linearize at the *current IPM iterate* (`kkt.current_x`), not at the
    # Krylov direction `xp` — the BM Jacobian is x-dependent (quadratic
    # constraint).
    NLPModels.hprod!(
        kkt.nlp, kkt.current_x, kkt.current_y, xp, kkt.hv_buf;
        obj_weight = one(eltype(yp)),
    )
    NLPModels.jtprod!(kkt.nlp, kkt.current_x, xd, kkt.jtv_buf)
    yp .= beta .* yp .+ alpha .* (kkt.hv_buf .+ kkt.jtv_buf)

    NLPModels.jprod!(kkt.nlp, kkt.current_x, xp, kkt.jv_buf)
    yd .= beta .* yd .+ alpha .* kkt.jv_buf
    return
end

# Variant called by MadNLP on `AbstractKKTVector` (extra `_kktmul!` at the
# end handles the bound-mult blocks; empty in our case but kept for safety).
function MadNLP.mul!(
    y::MadNLP.AbstractKKTVector,
    kkt::BMKKTSystem,
    x::MadNLP.AbstractKKTVector,
    alpha::Number, beta::Number,
)
    n, m = kkt.n, kkt.m
    _x = MadNLP.full(x)
    _y = MadNLP.full(y)
    _kkt_apply!(
        view(_y, 1:n), view(_y, (n+1):(n+m)),
        kkt,
        view(_x, 1:n), view(_x, (n+1):(n+m)),
        alpha, beta,
    )
    MadNLP._kktmul!(
        y, x,
        kkt.reg, kkt.du_diag, kkt.l_lower, kkt.u_lower, kkt.l_diag, kkt.u_diag,
        alpha, beta,
    )
    return y
end

# Variant called by `Krylov.kmul!` on plain `Vector{T}` (length `n + m`).
# Must mirror the `AbstractKKTVector` variant — including the `reg`/`du_diag`
# contributions that the sparse path applies via `_kktmul!`. Without these,
# MadNLP's `del_w` regularization never reaches Krylov, the curvature test
# can never be satisfied, and the IPM bails into restoration.
function LinearAlgebra.mul!(
    y::AbstractVector,
    kkt::BMKKTSystem,
    x::AbstractVector,
    alpha::Number, beta::Number,
)
    n, m = kkt.n, kkt.m
    yp = view(y, 1:n);       yd = view(y, (n+1):(n+m))
    xp = view(x, 1:n);       xd = view(x, (n+1):(n+m))
    _kkt_apply!(yp, yd, kkt, xp, xd, alpha, beta)
    # Match `_kktmul!` from the `AbstractKKTVector` path: add `reg` to the
    # primal block and `du_diag` to the dual block. These reflect MadNLP's
    # cumulative `del_w`/`del_c` Hessian/Jacobian perturbations.
    yp .+= alpha .* kkt.reg .* xp
    yd .+= alpha .* kkt.du_diag .* xd
    return y
end

LinearAlgebra.mul!(y::AbstractVector, kkt::BMKKTSystem, x::AbstractVector) =
    LinearAlgebra.mul!(y, kkt, x, true, false)

# Hessian-block matvec used by `InertiaFree`'s curvature test (`curv_test`
# in `MadNLP/src/IPM/solver.jl:785`). Computes `wx = (H + pr_diag) · t`
# where `t`, `wx` ∈ ℝⁿ are *primal only* (see `build_inertia_corrector`).
function MadNLP.mul_hess_blk!(
    wx::AbstractVector, kkt::BMKKTSystem, t::AbstractVector,
)
    NLPModels.hprod!(
        kkt.nlp, kkt.current_x, kkt.current_y, t, kkt.hv_buf;
        obj_weight = one(eltype(wx)),
    )
    copyto!(wx, kkt.hv_buf)
    wx .+= t .* kkt.pr_diag
    return wx
end

# `Aᵀ x`, linearized at the current IPM iterate (stashed in `kkt.current_x`
# by our `eval_jac_wrapper!`). Until that fix, this passed `y` as the
# evaluation point of the Jacobian — fine if `c` is linear, but `c(X) =
# jprod(model, X, X)` for our BM model is *quadratic* in `X`, so J depends
# on the linearization point and the wrong one produced bogus duals once
# `y` had non-trivial magnitude (feasibility happens to keep `y ≈ 0` so the
# bug stayed dormant; `max γ` revealed it).
function MadNLP.jtprod!(y::AbstractVector, kkt::BMKKTSystem, x::AbstractVector)
    NLPModels.jtprod!(kkt.nlp, kkt.current_x, x, y)
    return y
end

import SolverCore

# `BurerMonteiro.Solver`'s `solve!` calls
# `SolverCore.solve!(madnlp_solver, bm_model, stats; kws...)`. MadNLP's own
# entrypoints are `solve!(solver)` or `solve!(nlp, solver, stats)`; the
# (solver, nlp, stats) ordering doesn't exist. Bridge here.
#
# CompressedSensingIPM-style: only run `MadNLP.regular!`, never `robust!`.
# That bypasses restoration entirely — which is fine for us because
# restoration would call `eval_lag_hess_wrapper!`/`eval_jac_wrapper!` (which
# we stub to no-ops) and bail anyway.
# MadNLP `Status` → SolverCore status symbol (consumed by NLPModelsJuMP's
# `TERMINATION_STATUS`; `:unknown` → `MOI.OPTIMIZE_NOT_CALLED`, which makes
# JuMP throw `OptimizeNotCalled` on `value(...)` — avoid).
function _madnlp_to_solvercore_status(s::MadNLP.Status)
    s == MadNLP.SOLVE_SUCCEEDED                    && return :first_order
    s == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL         && return :acceptable
    s == MadNLP.SEARCH_DIRECTION_BECOMES_TOO_SMALL && return :small_step
    s == MadNLP.DIVERGING_ITERATES                 && return :unbounded
    s == MadNLP.INFEASIBLE_PROBLEM_DETECTED        && return :infeasible
    s == MadNLP.MAXIMUM_ITERATIONS_EXCEEDED        && return :max_iter
    s == MadNLP.MAXIMUM_WALLTIME_EXCEEDED          && return :max_time
    s == MadNLP.USER_REQUESTED_STOP                && return :user
    s == MadNLP.ERROR_IN_STEP_COMPUTATION          && return :neg_pred
    # In-progress states leaking out (e.g. we returned mid-iteration after the
    # line search failed) — treat as a stalled / slow-progress termination so
    # JuMP still sees a result and can read the iterate.
    s in (MadNLP.INITIAL, MadNLP.REGULAR, MadNLP.RESTORE, MadNLP.ROBUST,
          MadNLP.LINESEARCH_SUCCEEDED) && return :stalled
    return :exception
end

function SolverCore.solve!(
    solver::MadNLP.MadNLPSolver,
    nlp::NLPModels.AbstractNLPModel,
    stats::SolverCore.GenericExecutionStats;
    kws...,
)
    @info "Custom SolverCore.solve! → regular! only (skipping restoration)"
    MadNLP.print_init(solver)
    MadNLP.initialize!(solver)
    try
        MadNLP.regular!(solver)
    catch err
        @info "MadNLP.regular! threw — accepting current state" err
    end
    res = MadNLP.MadNLPExecutionStats(solver)
    stats.solution .= res.solution
    stats.objective = res.objective
    stats.iter      = res.iter
    SolverCore.set_status!(stats, _madnlp_to_solvercore_status(res.status))
    stats.dual_feas = res.dual_feas
    stats.primal_feas = res.primal_feas
    stats.solver_specific[:madnlp] = res
    return stats
end

# Probe `diag(H)` via `n` Hessian-vector products on standard basis vectors.
# Cheap for our smoke-test (`n = 14`); for the trigonometric `d = 100`
# benchmark this is `n ≈ 800` FFTs per IPM iter, still negligible compared
# to the Krylov inner loop.
function _refresh_preconditioner!(kkt::BMKKTSystem{QLP,T}) where {QLP,T}
    n, m = kkt.n, kkt.m
    e = kkt.precond_probe
    fill!(e, zero(T))
    for i in 1:n
        e[i] = one(T)
        NLPModels.hprod!(
            kkt.nlp, kkt.current_x, kkt.current_y, e, kkt.hv_buf;
            obj_weight = one(T),
        )
        # `|diag(H)| + reg + pr_diag`, guarded away from zero. Taking abs
        # keeps the preconditioner SPD even where H has indefinite diagonal
        # entries — for MINRES-QLP a SPD preconditioner is required.
        kkt.precond_diag[i] = max(abs(kkt.hv_buf[i]) + kkt.reg[i] + kkt.pr_diag[i], 1e-12)
        e[i] = zero(T)
    end
    for j in 1:m
        kkt.precond_diag[n + j] = max(abs(kkt.du_diag[j]), one(T))
    end
    return
end

function MadNLP.solve_kkt!(kkt::BMKKTSystem, w::MadNLP.AbstractKKTVector)
    MadNLP.reduce_rhs!(kkt, w)
    # The `[primal; dual]` block — what MadNLP's standard
    # `solve!(::AbstractReducedKKTSystem, w)` hands to its linear solver.
    b = MadNLP.primal_dual(w)
    rhs_norm = LinearAlgebra.norm(b)
    rhs_copy = copy(b)
    _refresh_preconditioner!(kkt)
    M = BMJacobiPrecond(kkt.precond_diag)
    Krylov.krylov_solve!(
        kkt.linear_solver, kkt, b;
        M = M, ldiv = true,
        atol = 1e-10, rtol = 1e-8, itmax = 10 * (kkt.n + kkt.m),
        verbose = 0,
    )
    x = Krylov.solution(kkt.linear_solver)
    # Compute true residual ‖K·x − b‖ against the operator we just gave Krylov
    Kx = similar(x)
    LinearAlgebra.mul!(Kx, kkt, x, true, false)   # routes through `_kkt_apply!`
    true_res = LinearAlgebra.norm(Kx .- rhs_copy)
    nit = Krylov.iteration_count(kkt.linear_solver)
    issolved = Krylov.issolved(kkt.linear_solver)
    kstatus = Krylov.statistics(kkt.linear_solver).status
    println(stderr,
        "Krylov solve: rhs=", rhs_norm,
        "  sol=", LinearAlgebra.norm(x),
        "  res=", true_res,
        "  it=", nit,
        "  solved=", issolved,
        "  ", kstatus,
    )
    copyto!(b, x)
    push!(kkt.krylov_iterations, nit)
    push!(kkt.krylov_residuals, Krylov.elapsed_time(kkt.linear_solver))
    MadNLP.finish_aug_solve!(kkt, w)
    return w
end
