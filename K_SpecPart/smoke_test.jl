# Pure-Julia regression smoke test for K-SpecPart.
#
# Exercises the internal pipeline that does NOT require the external binaries
# (gpmetis / hmetis / ilp_part / openroad):
#   read -> isolate_islands -> hypergraph2graph -> solve_eigs ->
#   reweigh_graph -> construct_tree -> distill_cuts_on_tree ->
#   two_way_linear_tree_sweep -> golden_evaluator -> overlay
#
# It prints a compact, deterministic summary so the same script can be run
# before and after refactoring to confirm parity.
#
# Usage:
#   julia smoke_test.jl [path/to/hypergraph.hgr]

include("specpart.jl")
using .SpecPart
using Random

const SP = SpecPart

function run_smoke(hgr_file::String; seed::Int = 0, ub::Int = 2,
                   nev::Int = 2, solver_iters::Int = 20, cycles::Int = 1)
    Random.seed!(seed)
    println("="^60)
    println("SMOKE TEST: ", hgr_file)
    println("="^60)

    hg = SP.read_hypergraph_file(hgr_file)
    println("read: vertices=", hg.num_vertices, " hyperedges=", hg.num_hyperedges)

    (phg, orig_idx, new_idx, unused_idx) = SP.isolate_islands(hg)
    println("isolated: vertices=", phg.num_vertices, " hyperedges=", phg.num_hyperedges)

    # Deterministic balanced 2-way hint: first half vs. second half.
    half = cld(phg.num_vertices, 2)
    hint = vcat(zeros(Int, half), ones(Int, phg.num_vertices - half))
    side0 = findall(==(0), hint)
    side1 = findall(==(1), hint)
    (hint_cut, hint_bal) = SP.golden_evaluator(phg, 2, hint)
    println("hint: cut=", hint_cut, " balance=", hint_bal)

    adj = SP.hypergraph2graph(phg, cycles)
    println("graphified: nnz=", length(adj.nzval))

    pindex = SP.__pindex__(side0, side1)
    X = SP.solve_eigs(phg, adj, pindex, false, nev, solver_iters)
    println("embedding: size=", size(X), " finite=", all(isfinite, X))

    # First eigenvector only, MST tree (type 2): pure-Julia, deterministic enough.
    Xv = X[:, 1:1]
    g = SP.reweigh_graph(adj, Xv, false)
    (tree, tree_matrix) = SP.construct_tree(g, Xv, 2)
    println("tree: nv=", SP.SimpleWeightedGraphs.nv(tree),
            " ne=", SP.SimpleWeightedGraphs.ne(tree))

    empty_p = SP.__pindex__(Int[], Int[])
    dc = SP.distill_cuts_on_tree(phg, empty_p, tree)

    total = sum(phg.vwts)
    # Relaxed capacities so the sweep reliably exercises its success path
    # (connected-components -> labels) and yields a deterministic cut number.
    maxcap = Int(ceil(total * (50 + 45) / 100))
    caps = [total - maxcap, maxcap]
    (part, sweep_cut, cutpt) = SP.two_way_linear_tree_sweep(tree, dc, phg, caps, 2, 2)
    println("tree_sweep: cost=", sweep_cut, " cutpoint=", cutpt)
    eval_cut = -1
    if cutpt > -1 && all(>=(0), part)
        (eval_cut, eval_bal) = SP.golden_evaluator(phg, 2, part)
        println("tree_sweep eval: cut=", eval_cut, " balance=", eval_bal)
    else
        part = copy(hint)
        println("tree_sweep eval: no valid cut (using hint for overlay)")
    end

    # Overlay the hint and the (tree-sweep or fallback) partition.
    (chg, clusters) = SP.overlay([hint, part], phg)
    println("overlay: clusters=", chg.num_vertices, " hyperedges=", chg.num_hyperedges)

    println("="^60)
    println("SUMMARY hint_cut=", hint_cut,
            " sweep_cut=", eval_cut,
            " overlay_clusters=", chg.num_vertices)
    println("="^60)
    return (hint_cut, eval_cut, chg.num_vertices)
end

hgr = length(ARGS) >= 1 ? ARGS[1] : "ilp_partitioner/test/ibm10.hgr"
run_smoke(hgr)
