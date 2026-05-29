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
end

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
    copyto!(kkt.current_y, l)
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
function LinearAlgebra.mul!(
    y::AbstractVector,
    kkt::BMKKTSystem,
    x::AbstractVector,
    alpha::Number, beta::Number,
)
    n, m = kkt.n, kkt.m
    _kkt_apply!(
        view(y, 1:n), view(y, (n+1):(n+m)),
        kkt,
        view(x, 1:n), view(x, (n+1):(n+m)),
        alpha, beta,
    )
    return y
end

LinearAlgebra.mul!(y::AbstractVector, kkt::BMKKTSystem, x::AbstractVector) =
    LinearAlgebra.mul!(y, kkt, x, true, false)

# `Aᵀ x`, just routes to the NLP's `jtprod!`.
function MadNLP.jtprod!(y::AbstractVector, kkt::BMKKTSystem, x::AbstractVector)
    NLPModels.jtprod!(kkt.nlp, y, x, y)
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

function MadNLP.solve_kkt!(kkt::BMKKTSystem, w::MadNLP.AbstractKKTVector)
    MadNLP.reduce_rhs!(kkt, w)
    # The `[primal; dual]` block — what MadNLP's standard
    # `solve!(::AbstractReducedKKTSystem, w)` hands to its linear solver.
    b = MadNLP.primal_dual(w)
    rhs_norm = LinearAlgebra.norm(b)
    rhs_copy = copy(b)
    Krylov.krylov_solve!(
        kkt.linear_solver, kkt, b;
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
