# Required branches
#
#   SumOfSquares       → branch `bl/fft`    (LowRankOpt + FFT wiring)
#   LowRankOpt         → branch `bl/sampling` (incl. `hprod!` sign fix)
#   MultivariateBases  → branch `bl/trigpolys` (TrigEvalMatrix + batched mul!)
#
# Run from this directory:
#
#   julia --project=. bm_madnlp_bench.jl

using SumOfSquares
import DynamicPolynomials
import MultivariateBases as MB
import LowRankOpt as LRO
import MadNLP
import Krylov
import NLPModels
import Dualization
import LinearAlgebra
import MathOptInterface as MOI
import Random

include(joinpath(@__DIR__, "bm_madnlp_kkt.jl"))

DynamicPolynomials.@polyvar x

# Custom sub_solver factory that builds a MadNLPSolver with our BMKKTSystem.
function MadNLPMinresSolver(nlp::NLPModels.AbstractNLPModel; qlp::Bool = true, kws...)
    # Deterministic init so every run sees the same KKT trace.
    Random.seed!(0)
    # Override the random `meta.x0` with a feasibility-friendlier starting
    # point. Starting from `rand(n)/n` (small) makes the Hessian near-singular
    # at iter 0, which gives Newton steps that overshoot the optimum and the
    # IPM oscillates. Use `rand(n)` (BMSOSAL's choice) for a non-degenerate
    # initial Hessian.
    n = length(nlp.meta.x0)
    nlp.meta.x0 .= rand(eltype(nlp.meta.x0), n)
    # Quick check whether the BM objective is identically zero (which would
    # make the dual-Newton step purely feasibility-driven).
    f0 = NLPModels.obj(nlp, nlp.meta.x0)
    g0 = NLPModels.grad(nlp, nlp.meta.x0)
    @info "BM model at x0" obj_x0=f0 grad_norm=LinearAlgebra.norm(g0) nvar=n ncon=length(nlp.meta.y0)
    return MadNLP.MadNLPSolver(
        nlp;
        callback = MadNLP.DenseCallback,
        kkt_system = BMKKTSystem{qlp},
        linear_solver = MadNLP.LapackCPUSolver,  # unused — Krylov takes over
        nlp_scaling = false,                     # `jac_dense!` not implemented
        # We can't honestly report KKT inertia from a matrix-free MINRES
        # solve, so use MadNLP's inertia-free correction (Curtis–Schenk–Wächter
        # heuristic). It just needs two extra `solve_kkt!`s per Newton step
        # to validate that the computed direction is a descent step.
        inertia_correction_method = MadNLP.InertiaFree,
        # Stabilisation : la `curv_test` par défaut tolère ζ=0 (curvature ≥ 0),
        # ce qui laisse passer des directions presque-singulières → pas de
        # Newton qui dépassent l'optimum. Demander une curvature strictement
        # positive et un régulariseur primal de base force un comportement
        # type «Levenberg-Marquardt léger».
        inertia_free_tol = 1e-8,
        default_primal_regularization   = 1e-8,
        default_dual_regularization     = 1e-8,
        # Default least-squares dual init solves K·[d;y] = [-g;-c] once,
        # which on our singular KKT can throw `y` to wild values. Starting
        # from `y = 0` lets MadNLP build up the duals iteration by iteration.
        dual_initialization_method = MadNLP.DualInitializeSetZero,
        print_level = MadNLP.INFO,
        max_iter = 500,
    )
end

inner = optimizer_with_attributes(
    LRO.Optimizer,
    "solver"          => LRO.BurerMonteiro.Solver,
    "sub_solver"      => MadNLPMinresSolver,
    "ranks"           => [4],
    "square_scalars"  => true,
)
bmlbfgs = Dualization.dual_optimizer(inner; assume_min_if_feasibility = true)

function with_lro_bridges!(model)
    backend = JuMP.backend(model)
    SumOfSquares.Bridges.add_all_bridges(backend.optimizer, Float64)
    MOI.Bridges.remove_bridge(
        backend.optimizer,
        SumOfSquares.Bridges.Constraint.ImageBridge{Float64},
    )
end

p = x^4 - 4x^3 - 2x^2 + 12x + 3

function smoke(feas::Bool)
    println("\n== BMKKTSystem smoke: ", feas ? "feasibility (γ=-6)" : "max γ", " ==")
    model = Model(bmlbfgs)
    set_silent(model)
    γ = -6
    if !feas
        @variable(model, γ)
        @objective(model, Max, γ)
    end
    with_lro_bridges!(model)
    @constraint(model, p - γ in SOSCone(), zero_basis = MB.BoxSampling([-1.0], [1.0]))
    optimize!(model)
    println("primal_status = ", primal_status(model))
    if !feas
        println("value(γ)      = ", value(γ), "    (expected ≈ -6)")
    end
end

smoke(true)
smoke(false)

# Cross-check : same problem, but with Percival as the sub_solver. If Percival
# *also* diverges on `max γ`, the bug is upstream (Dualization / bridges),
# not in our MadNLP integration.
println("\n== Percival cross-check (max γ) ==")
import Percival
percival_inner = optimizer_with_attributes(
    LRO.Optimizer,
    "solver"         => LRO.BurerMonteiro.Solver,
    "sub_solver"     => Percival.PercivalSolver,
    "ranks"          => [4],
    "square_scalars" => true,
)
percival_bmlbfgs = Dualization.dual_optimizer(percival_inner; assume_min_if_feasibility = true)
let
    model = Model(percival_bmlbfgs)
    set_silent(model)
    @variable(model, γ)
    @objective(model, Max, γ)
    with_lro_bridges!(model)
    @constraint(model, p - γ in SOSCone(), zero_basis = MB.BoxSampling([-1.0], [1.0]))
    optimize!(model)
    println("primal_status = ", primal_status(model))
    println("value(γ)      = ", value(γ), "    (expected ≈ -6)")
end

# Should give
# == BMKKTSystem smoke: max γ s.t. p − γ ∈ SOS ==
# ┌ Info: BM model at x0
# │   obj_x0 = 0.00140682890105754
# │   grad_norm = 0.12033747127060285
# │   nvar = 14
# └   ncon = 5
# [ Info: Custom SolverCore.solve! → regular! only (skipping restoration)
# Number of nonzeros in constraint Jacobian............:       70
# Number of nonzeros in Lagrangian Hessian.............:      105
# 
# Total number of variables............................:       14
#                      variables with only lower bounds:        0
#                 variables with lower and upper bounds:        0
#                      variables with only upper bounds:        0
# Total number of equality constraints.................:        5
# Total number of inequality constraints...............:        0
#         inequality constraints with only lower bounds:        0
#    inequality constraints with lower and upper bounds:        0
#         inequality constraints with only upper bounds:        0
# 
# iter    objective    inf_pr   inf_du inf_compl lg(mu) lg(rg) alpha_pr ir ls
#    0  1.4068289e-03 9.97e+00 1.00e-01 0.00e+00  -1.0     -   0.00e+00 10  0 
# primal_status = NO_SOLUTION
