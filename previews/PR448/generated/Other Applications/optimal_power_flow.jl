using SumOfSquares

import PowerModels
import PGLib
PowerModels.silence()

data = PGLib.pglib("pglib_opf_case14_ieee")

import Ipopt
local_solution = PowerModels.solve_ac_opf(
    data,
    optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
)
AC = local_solution["objective"]

using DynamicPolynomials

bfind(v, x) = findfirst(isequal(x), v)
fl_sum(itr) = mapreduce(identity, +, itr, init = 0.0)

function mkpoly(x, supp, coe)
    return sum(
        coe[t] * prod(j -> x[j], supp[t]; init = one(eltype(coe))) for
        t in eachindex(coe)
    )
end
normalize_poly(p) = p / maximum(abs, coefficients(p))

function build_opf_pop(data; AngleCons = true, LineLimit = "relax")
    PowerModels.standardize_cost_terms!(data, order = 2)
    ref = PowerModels.build_ref(data)[:it][PowerModels.pm_it_sym][:nw][0]
    nbus = length(ref[:bus])
    ng = length(ref[:gen])
    n = 2 * nbus + 2 * ng
    @polyvar x[1:n]
    ineqs = Any[]
    eqs = Any[]
    gens = sort!(collect(keys(ref[:gen])))
    bus = sort!(collect(keys(ref[:bus])))

    # Objective: the sum of the (quadratic) generation costs.
    coe = Float64[sum(gen["cost"][3] for (_, gen) in ref[:gen])]
    supp = Vector{Int}[Int[]]
    for i in 1:ng
        gen = ref[:gen][gens[i]]
        push!(coe, gen["cost"][2], gen["cost"][1])
        push!(supp, [2nbus + i], [2nbus + i, 2nbus + i])
    end
    f = mkpoly(x, supp, coe)

    # Voltage magnitude: vmin² ≤ eᵢ² + fᵢ² ≤ vmax².
    for i in 1:nbus
        push!(ineqs, mkpoly(x, Vector{Int}[[], [i, i], [i+nbus, i+nbus]],
            [-ref[:bus][bus[i]]["vmin"]^2, 1, 1]))
        push!(ineqs, mkpoly(x, Vector{Int}[[], [i, i], [i+nbus, i+nbus]],
            [ref[:bus][bus[i]]["vmax"]^2, -1, -1]))
    end

    # Angle differences and (relaxed) thermal limits, one pair per branch.
    for (_, branch) in ref[:branch]
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]; b_fr = branch["b_fr"]
        g_to = branch["g_to"]; b_to = branch["b_to"]
        tm = branch["tap"]
        vr = bfind(bus, branch["f_bus"])
        vt = bfind(bus, branch["t_bus"])
        srt = sort([vr, vt])
        if AngleCons
            p1 = mkpoly(x, Vector{Int}[srt, srt .+ nbus, [vt, vr+nbus], [vr, vt+nbus]],
                [tan(branch["angmax"]), tan(branch["angmax"]), -1, 1])
            p2 = mkpoly(x, Vector{Int}[[vt, vr+nbus], [vr, vt+nbus], srt, srt .+ nbus],
                [1, -1, -tan(branch["angmin"]), -tan(branch["angmin"])])
            push!(ineqs, normalize_poly(p1), normalize_poly(p2))
        end
        if LineLimit == "relax"
            ab1 = (g+g_fr)^2 + (b+b_fr)^2
            cd1 = (-g*tr+b*ti)^2 + (b*tr+g*ti)^2
            acbd1 = (g+g_fr)*(-g*tr+b*ti) - (b+b_fr)*(b*tr+g*ti)
            bcad1 = -(b+b_fr)*(-g*tr+b*ti) - (g+g_fr)*(b*tr+g*ti)
            ab2 = (g+g_to)^2*tm^4 + (b+b_to)^2*tm^4
            cd2 = (g*tr+b*ti)^2 + (-b*tr+g*ti)^2
            acbd2 = -(g+g_to)*tm^2*(g*tr+b*ti) + (b+b_to)*tm^2*(-b*tr+g*ti)
            bcad2 = (b+b_to)*tm^2*(g*tr+b*ti) + (g+g_to)*tm^2*(-b*tr+g*ti)
            mvr = ref[:bus][bus[vr]]["vmin"]^2
            mvt = ref[:bus][bus[vt]]["vmin"]^2
            p1 = mkpoly(x, Vector{Int}[[], [vr, vr], [vr+nbus, vr+nbus], [vt, vt], [vt+nbus, vt+nbus], srt, [vr, vt+nbus], [vt, vr+nbus], srt .+ nbus],
                [branch["rate_a"]^2*tm^4/mvr, -ab1, -ab1, -cd1, -cd1, -2acbd1, 2bcad1, -2bcad1, -2acbd1])
            p2 = mkpoly(x, Vector{Int}[[], [vt, vt], [vt+nbus, vt+nbus], [vr, vr], [vr+nbus, vr+nbus], srt, [vt, vr+nbus], [vr, vt+nbus], srt .+ nbus],
                [branch["rate_a"]^2*tm^4/mvt, -ab2, -ab2, -cd2, -cd2, -2acbd2, 2bcad2, -2bcad2, -2acbd2])
            push!(ineqs, normalize_poly(p1), normalize_poly(p2))
        end
    end

    # Generation bounds: pmin ≤ pᵍ ≤ pmax and qmin ≤ qᵍ ≤ qmax.
    for i in 1:ng
        gen = ref[:gen][gens[i]]
        p = mkpoly(x, Vector{Int}[[], [2nbus+i], [2nbus+i, 2nbus+i]],
            [-gen["pmin"]*gen["pmax"], gen["pmin"]+gen["pmax"], -1])
        push!(ineqs, normalize_poly(p))
        p = mkpoly(x, Vector{Int}[[], [2nbus+ng+i], [2nbus+ng+i, 2nbus+ng+i]],
            [-gen["qmin"]*gen["qmax"], gen["qmin"]+gen["qmax"], -1])
        push!(ineqs, normalize_poly(p))
    end

    # Power-flow balance at each bus (Kirchhoff's laws): equality constraints.
    for (r, i) in enumerate(bus)
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        na = 4 * length(ref[:bus_arcs][i]) + 3
        coe1 = zeros(na); supp1 = Vector{Vector{Int}}(undef, na)
        coe2 = zeros(na); supp2 = Vector{Vector{Int}}(undef, na)
        supp1[1:3] = [Int[], [r, r], [r+nbus, r+nbus]]
        supp2[1:3] = [Int[], [r, r], [r+nbus, r+nbus]]
        coe1[1] = fl_sum(load["pd"] for load in bus_loads)
        coe2[1] = fl_sum(load["qd"] for load in bus_loads)
        sgs = fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs = fl_sum(shunt["bs"] for shunt in bus_shunts)
        coe1[2:3] = [sgs, sgs]
        coe2[2:3] = [-sbs, -sbs]
        j = 1
        for flow in ref[:bus_arcs][i]
            branch = ref[:branch][flow[1]]
            vr = bfind(bus, branch["f_bus"])
            vt = bfind(bus, branch["t_bus"])
            srt = sort([vr, vt])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]; b_fr = branch["b_fr"]
            g_to = branch["g_to"]; b_to = branch["b_to"]
            tm = branch["tap"]
            a1 = (g+g_fr)/tm^2; b1 = -(b+b_fr)/tm^2
            c1 = (-g*tr+b*ti)/tm^2; d1 = (b*tr+g*ti)/tm^2
            a2 = g+g_to; b2 = -(b+b_to)
            c2 = -(g*tr+b*ti)/tm^2; d2 = -(-b*tr+g*ti)/tm^2
            supp1[j+3:j+6] = [srt, [vt, vr+nbus], [vr, vt+nbus], srt .+ nbus]
            supp2[j+3:j+6] = [srt, [vt, vr+nbus], [vr, vt+nbus], srt .+ nbus]
            if vr == r
                coe1[2:3] .+= a1; coe1[j+3:j+6] = [c1, -d1, d1, c1]
                coe2[2:3] .+= b1; coe2[j+3:j+6] = [d1, c1, -c1, d1]
            else
                coe1[2:3] .+= a2; coe1[j+3:j+6] = [c2, d2, -d2, c2]
                coe2[2:3] .+= b2; coe2[j+3:j+6] = [d2, -c2, c2, d2]
            end
            j += 4
        end
        for gen_id in ref[:bus_gens][i]
            gen = bfind(gens, gen_id)
            push!(supp1, [2nbus + gen]); push!(coe1, -1)
            push!(supp2, [2nbus + ng + gen]); push!(coe2, -1)
        end
        push!(eqs, mkpoly(x, supp1, coe1))
        push!(eqs, mkpoly(x, supp2, coe2))
    end

    # Reference bus: the imaginary part of the voltage is set to zero.
    for key in keys(ref[:ref_buses])
        i = bfind(bus, key)
        push!(eqs, mkpoly(x, Vector{Int}[[i+nbus, i+nbus]], Float64[1]))
    end
    return x, f, ineqs, eqs
end

x, f, ineqs, eqs = build_opf_pop(data)
length(x), length(ineqs), length(eqs)

import Clarabel

function build_domain(ineqs, eqs)
    gs = vcat(ineqs, eqs, [-h for h in eqs])
    return mapreduce(g -> (@set g >= 0), intersect, gs)
end

function relaxation(f, ineqs, eqs; sparsity = Sparsity.NoPattern())
    model = SOSModel(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, t)
    @objective(model, Max, t)
    con_ref = @constraint(
        model, f >= t,
        domain = build_domain(ineqs, eqs), maxdegree = 2, sparsity = sparsity,
    )
    optimize!(model)
    return model, con_ref
end

function block_sizes(con_ref)
    g = gram_matrix(con_ref)
    if g isa SumOfSquares.BlockDiagonalGramMatrix
        return sort!([length(b.basis) for b in g.blocks]; rev = true)
    else
        return [length(g.basis)]
    end
end

model, con_ref = relaxation(f, ineqs, eqs)
dense_bound = objective_value(model)
solution_summary(model)

dense_bound

gap(opt) = 100 * (AC - opt) / AC
gap(dense_bound)

dense_blocks = block_sizes(con_ref)
maximum(dense_blocks)

model, con_ref = relaxation(f, ineqs, eqs; sparsity = Sparsity.Variable())
correlative_bound = objective_value(model)
gap(correlative_bound)

correlative_blocks = block_sizes(con_ref)
correlative_blocks

model, con_ref = relaxation(f, ineqs, eqs; sparsity = Sparsity.Monomial(ChordalCompletion()))
term_bound = objective_value(model)
gap(term_bound)

term_blocks = block_sizes(con_ref)
term_blocks

using Printf
for (name, bound, blocks) in [
    ("dense", dense_bound, dense_blocks),
    ("correlative", correlative_bound, correlative_blocks),
    ("term", term_bound, term_blocks),
]
    @printf(
        "%-12s bound = %.2f   gap = %6.3f%%   #blocks = %3d   max block = %d\n",
        name, bound, gap(bound), length(blocks), maximum(blocks),
    )
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
