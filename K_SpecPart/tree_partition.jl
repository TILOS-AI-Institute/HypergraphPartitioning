using Combinatorics
using Laplacians 
using SimpleGraphs
using Clustering
include("cut_distillation.jl")
include("metis.jl")
include("degree_aware_prims.jl")
include("extract_hypergraph.jl")

# Number of k-means restarts used to generate direct (non-tree) spectral
# clustering candidates from the embedding. Default 0 (opt-in via
# KSPECPART_KMEANS_SEEDS): on the tested designs the balance-repaired k-means
# candidates had much higher cuts than the tree candidates and did not improve
# results, while adding FM/runtime cost. Code retained for experimentation.
const NUM_KMEANS_SEEDS = parse(Int, get(ENV, "KSPECPART_KMEANS_SEEDS", "0"))

# Greedy capacity-respecting assignment of points (rows of distance matrix D,
# n x num_parts, smaller = closer) to blocks so no block exceeds `maxcap`.
# Vertices that are most decisive (smallest min distance) are placed first.
function balanced_assign(D::Matrix{Float64}, vwts::Vector{Int},
                         num_parts::Int, maxcap::Int)
    n = size(D, 1)
    part = fill(-1, n)
    load = zeros(Int, num_parts)
    order = sortperm([minimum(@view D[i, :]) for i in 1:n])
    for i in order
        prefs = sortperm(@view D[i, :])
        placed = false
        for b in prefs
            if load[b] + vwts[i] <= maxcap
                part[i] = b - 1
                load[b] += vwts[i]
                placed = true
                break
            end
        end
        if !placed
            b = argmin(load)
            part[i] = b - 1
            load[b] += vwts[i]
        end
    end
    return part
end

# Direct spectral-clustering candidates: cluster the embedding into num_parts
# groups (k-means), repair to balance, and return `[partition, cutsize]` entries.
# These are a distinct candidate class from the tree sweeps; they are FM-refined
# downstream like any other candidate.
function kmeans_candidates(X::Array{Float64}, hgraph::__hypergraph__,
                           num_parts::Int, capacities::Vector{Int}, seed::Int)
    cands = []
    NUM_KMEANS_SEEDS <= 0 && return cands
    n = size(X, 1)
    d = size(X, 2)
    (n <= num_parts || num_parts < 2) && return cands
    maxcap = capacities[2]
    Xt = permutedims(X)          # d x n, as Clustering.kmeans expects
    for s in 1:NUM_KMEANS_SEEDS
        Random.seed!(hash((seed, :kmeans, s)))
        local result
        try
            result = kmeans(Xt, num_parts; maxiter=100)
        catch
            continue
        end
        centers = result.centers  # d x num_parts
        D = zeros(Float64, n, num_parts)
        @inbounds for b in 1:num_parts
            for i in 1:n
                acc = 0.0
                for q in 1:d
                    diff = X[i, q] - centers[q, b]
                    acc += diff * diff
                end
                D[i, b] = acc
            end
        end
        part = balanced_assign(D, hgraph.vwts, num_parts, maxcap)
        (cut, ~) = golden_evaluator(hgraph, num_parts, part)
        @debug "k-means candidate (seed $s) cut=$cut"
        push!(cands, [part, cut])
    end
    return cands
end

function reweigh_graph(adj::SparseMatrixCSC, 
                    X::AbstractArray, 
                    lst::Bool)
    n = size(X, 1)
    nev = size(X, 2)
    offset = 0
    g = SimpleWeightedGraph(adj)
    g_copy = deepcopy(g)
    ewts = g_copy.weights
    for i in 1:length(ewts.colptr)-1
        row_len = ewts.colptr[i+1] - ewts.colptr[i]
        for row in 1:row_len
            row_x = ewts.rowval[row+offset]
            distance = 0.0
            for d in 1:nev
                span = X[row_x, d] - X[i, d]
                if lst == true
                    if span == 0.0
                        distance += 1e09
                    else
                        distance += 1/(span*span) 
                    end
                else
                    distance += abs(span)
                end
            end
            ewts.nzval[row+offset] = distance
        end
        offset += row_len
    end
    g_copy.weights = ewts
    return g_copy
end

function reweigh_graph_with_cuts(adj::SparseMatrixCSC, 
                                hgraph::__hypergraph__,
                                X::AbstractArray, 
                                lst::Bool)
    n = size(X, 1)
    nev = size(X, 2)
    offset = 0
    g = SimpleWeightedGraph(adj)
    g_copy = deepcopy(g)
    ewts = g_copy.weights
    (~, ~, he_wts) = findnz(ewts)
    (min_wt, max_wt) = extrema(he_wts)
    norm_factor = max_wt - min_wt
    norm_factor == 0 && (norm_factor = 1.0)
    he_wts ./= norm_factor
    for i in 1:length(ewts.colptr)-1
        row_len = ewts.colptr[i+1] - ewts.colptr[i]
        for row in 1:row_len
            row_x = ewts.rowval[row+offset]
            distance = 0.0
            for d in 1:nev
                span = X[row_x, d] - X[i, d]
                if lst == true
                    if span == 0.0
                        distance += 1e09
                    else
                        distance += 1/(span*span) 
                    end
                else
                    distance += abs(span)
                end
            end
            ewts.nzval[row+offset] = distance
        end
        offset += row_len
    end
    (ii, jj, algebraic_wts) = findnz(ewts)
    (min_wt, max_wt) = extrema(algebraic_wts)
    norm_factor = max_wt - min_wt
    norm_factor == 0 && (norm_factor = 1.0)
    algebraic_wts ./= norm_factor
    he_factor = 100.0 
    algebraic_factor = 1000.0
    total_wts = he_factor .* he_wts + algebraic_factor .* algebraic_wts
    g_matrix = sparse(ii, jj, total_wts)
    return SimpleWeightedGraph(g_matrix)
end

function construct_tree(g::SimpleWeightedGraph, 
                        X::AbstractArray, 
                        tree_type::Int)
    tree = SimpleWeightedGraphs.SimpleGraph(g)
    n = SimpleWeightedGraphs.nv(g)
    tree_matrix = spzeros(n, n)
    if tree_type == 1   #lsst construction
        vtx_ids = Vector{Int}(1:n)
        fiedler = X[:,1]
        sorted_vtx_ids = sortperm(fiedler)
        reverse_map = zeros(Int, n)
        for i in 1:n
            reverse_map[sorted_vtx_ids[i]] = i
        end
        (i, j, w) = findnz(g.weights)
        tree_reordered_mat = sparse(sorted_vtx_ids[i], sorted_vtx_ids[j], w)
        tree_reordered_graph = SimpleWeightedGraphs.SimpleWeightedGraph(tree_reordered_mat)
        lsst = akpw(tree_reordered_graph.weights)
        (i, j, w) = findnz(lsst)
        lsst = sparse(reverse_map[i], reverse_map[j], w)
        tree = SimpleWeightedGraphs.SimpleGraph(lsst)
        tree_matrix = lsst
    elseif tree_type == 2   #mst construction
        #mst = SimpleWeightedGraphs.prim_mst(g)
        mst = degrees_aware_prim_mst(g, 10)
        i = zeros(Int, length(mst))
        j = zeros(Int, length(mst))
        w = zeros(length(mst))
        for k in 1:length(mst)
            vpair = mst[k]
            vsrc = vpair.src
            vdst = vpair.dst
            i[k] = vsrc
            j[k] = vdst
            w[k] = g.weights[vsrc, vdst]
        end
        tree = SimpleWeightedGraphs.SimpleGraph(mst)
        tree_matrix = sparse(i, j, w, n, n)
        tree_matrix += tree_matrix'
    elseif tree_type == 3   #path construction
        vtxs = sortperm(X)
        i = vtxs[1:end-1]
        j = vtxs[2:end]
        w = abs.(X[j] - X[i])
        w_z = findall(w .== 0.0)
        w[w_z] .= 1e-6
        tree_matrix = sparse(i, j, w, n, n)
        tree_matrix += tree_matrix'
        tree = SimpleWeightedGraphs.SimpleGraph(tree_matrix)
    else
        @warn "Please select correct tree type!"
    end
    return tree, tree_matrix
end

function two_way_linear_tree_sweep(T::SimpleWeightedGraphs.SimpleGraph, 
                                distilled_cuts::__cut_profile__,
                                hgraph::__hypergraph__, 
                                capacities::Vector{Int}, 
                                num_parts::Int,
                                solns::Int)
    hwts = hgraph.hwts
    num_vertices = hgraph.num_vertices
    num_hyperedges = hgraph.num_hyperedges
    vwts = hgraph.vwts
    partition_matrix = zeros(Int, solns, num_vertices)
    vtx_cuts = distilled_cuts.vtx_cuts
    edge_cuts = distilled_cuts.edge_cuts
    edge_diff = distilled_cuts.edge_diff
    pred = distilled_cuts.pred
    edge_terminators = distilled_cuts.edge_terminators
    p = distilled_cuts.p
    forced_type = distilled_cuts.forced_type
    forced_0 = distilled_cuts.forced_0
    forced_1 = distilled_cuts.forced_1
    forced_01 = distilled_cuts.forced_01
    FB0 = distilled_cuts.FB0
    FB1 = distilled_cuts.FB1
    edge_cuts_0 = distilled_cuts.edge_cuts_0
    edge_cuts_1 = distilled_cuts.edge_cuts_1
    nforced_0 = sum(hwts[forced_0])
    nforced_1 = sum(hwts[forced_1])
    nforced_01 = sum(hwts[forced_01])
    twt = sum(vwts)
    edge_cuts[1] = num_hyperedges
    vtx_cuts[1] = 0
    exc_0 = zeros(Int, num_vertices)
    exc_1 = zeros(Int, num_vertices)
    area_part = zeros(Int, 2)
    cut_cost_0 = zeros(num_vertices)
    cut_cost_1 = zeros(num_vertices)
    ratio_cost = zeros(num_vertices)
    cut_cost = zeros(num_vertices)
    area_cost = zeros(num_vertices)
    total_cost = zeros(num_vertices)
    polarity = zeros(Int, num_vertices)
    status_flag = zeros(Int, num_vertices)
    area_util_0 = zeros(Int, num_vertices)
    area_util_1 = zeros(Int, num_vertices)

    for i in 1:num_vertices
        exc_0[i] = edge_cuts[i] + nforced_0 - FB0[i] + edge_cuts_1[i] + nforced_01
        exc_1[i] = edge_cuts[i] + nforced_1 - FB1[i] + edge_cuts_0[i] + nforced_01
    end

    exc_0[1] = exc_1[1] = num_hyperedges

    for i in 1:num_vertices
        cut_cost_0[i] = exc_0[i]
        cut_cost_1[i] = exc_1[i]
        (cut_cost[i], pol) = findmin([cut_cost_0[i], cut_cost_1[i]])
        polarity[i] = pol-1

        if pol == 0
            area_util_0[i] = vtx_cuts[i]
            area_util_1[i] = twt - vtx_cuts[i]
        else
            area_util_1[i] = vtx_cuts[i]
            area_util_0[i] = twt - vtx_cuts[i]
        end

        #=SimpleWeightedGraphs.rem_edge!(T, i, pred[i])
        comps = SimpleWeightedGraphs.connected_components(T)
        cc = find_labels(comps, num_vertices)
        (cutsize, bal) = golden_evaluator(hgraph, 2, cc)
        println("[debug] bal comparison ", bal, ", ", area_util_0[i], " and ", area_util_1[i])
        SimpleWeightedGraphs.add_edge!(T, i, pred[i])=#

        if area_util_0[i] > capacities[2] || area_util_1[i] > capacities[2]
            area_cost[i] = 1e09
            #=area_cost_0 = area_util_0[i] > capacities[2] ? area_util_0[i] - capacities[2] : 0
            area_cost_1 = area_util_1[i] > capacities[2] ? area_util_1[i] - capacities[2] : 0
            area_cost[i] = (area_cost_0 + area_cost_1) * 1000
            if (i < 10) 
                println("Iteration ", i, " cutcost ", cut_cost[i], " areacost ", area_cost[i])
            end=#
        end

        ratio_cost[i] = cut_cost[i]/(area_util_0[i] * area_util_1[i])
        total_cost[i] = cut_cost[i] + area_cost[i]
    end

    cut_idxs = sortperm(total_cost)
    cut_point = cut_idxs[1]
    
    if total_cost[cut_idxs[1]] >= 1e09
        cut_idxs = sortperm(ratio_cost)
        cut_point = cut_idxs[1]
        overflow = false
        if area_util_0[cut_point] > capacities[2] || 
            area_util_1[cut_point] > capacities[2]
            overflow = true
        end
        i = 2
        while overflow == true
            cut_point = cut_idxs[i]
            overflow = false
            if area_util_0[cut_point] > capacities[2] || 
                area_util_1[cut_point] > capacities[2]
                overflow = true
            end
            i += 1
            if i > num_vertices
                cut_point= -1
                break
            end
        end
    end
    
    #println("[debug] cutpoint ", cut_point)
    partition = -ones(Int, hgraph.num_vertices)
    cutsize = 1e09
    if cut_point > -1
        SimpleWeightedGraphs.rem_edge!(T, cut_point, pred[cut_point])
        comps = SimpleWeightedGraphs.connected_components(T)
        partition = find_labels(comps, num_vertices)
        @debug "tree sweep cut=$(cut_cost[cut_point])"
        (cutsize, ~) = golden_evaluator(hgraph, num_parts, partition)
    else
        @debug "tree sweep found no valid balanced cut"
    end
    return (partition, cutsize, cut_point)
end

function METIS_tree_partition(T::SimpleWeightedGraphs.SimpleGraph, 
                            distilled_cuts::__cut_profile__,
                            hgraph::__hypergraph__, 
                            seed::Int,
                            metis_path::String,
                            metis_opts::Int,
                            num_parts::Int,
                            ub_factor::Real)
    hwts = hgraph.hwts
    num_vertices = hgraph.num_vertices
    num_hyperedges = hgraph.num_hyperedges
    vwts = hgraph.vwts
    vtx_cuts = distilled_cuts.vtx_cuts
    edge_cuts = distilled_cuts.edge_cuts
    edge_diff = distilled_cuts.edge_diff
    pred = distilled_cuts.pred
    edge_terminators = distilled_cuts.edge_terminators
    p = distilled_cuts.p
    forced_type = distilled_cuts.forced_type
    forced_0 = distilled_cuts.forced_0
    forced_1 = distilled_cuts.forced_1
    forced_01 = distilled_cuts.forced_01
    FB0 = distilled_cuts.FB0
    FB1 = distilled_cuts.FB1
    edge_cuts_0 = distilled_cuts.edge_cuts_0
    edge_cuts_1 = distilled_cuts.edge_cuts_1
    nforced_0 = sum(hwts[forced_0])
    nforced_1 = sum(hwts[forced_1])
    nforced_01 = sum(hwts[forced_01])
    exc_0 = zeros(Int, num_vertices)
    exc_1 = zeros(Int, num_vertices)
    cut_cost = zeros(num_vertices)
    T_matrix = sparse(T)

    for i in 1:hgraph.num_vertices
        if pred[i] == i
            cut_cost[i] = 1e09
            continue
        end
        exc_0[i] = edge_cuts[i] + nforced_0 - FB0[i] + edge_cuts_1[i] + nforced_01
        exc_1[i] = edge_cuts[i] + nforced_1 - FB1[i] + edge_cuts_0[i] + nforced_01
        cut_cost[i] = min(exc_0[i], exc_1[i])
        parent = pred[i]
        T_matrix[i, parent] = cut_cost[i]
        T_matrix[parent, i] = cut_cost[i]
    end
    g = SimpleWeightedGraph(T_matrix)
    gname = build_metis_graph(g, metis_opts)
    metis(metis_path, gname, num_parts, seed, ub_factor, metis_opts)
    pname = gname * ".part." * string(num_parts)
    partition = zeros(Int, hgraph.num_vertices)
    partition_i = 0
    for ln in eachline(pname)
        partition_i += 1
        partition[partition_i] = parse(Int, ln)
    end
    rm(pname, force=true)
    (cutsize, ~) = golden_evaluator(hgraph, num_parts, partition)
    @debug "METIS-on-tree cut=$cutsize"
    return (partition, cutsize)
end

function k_way_linear_tree_sweep(T::SimpleWeightedGraphs.SimpleGraph, 
                                distilled_cuts::__cut_profile__,
                                fixed_vertices::__pindex__,
                                hgraph::__hypergraph__, 
                                capacities::Vector{Int}, 
                                num_parts::Int,
                                solns::Int)
    recursion_levels = num_parts-1
    hypergraph_recursive = deepcopy(hgraph)
    tree_recursive = deepcopy(T)
    total_wts = sum(hgraph.vwts) 
    distilled_cuts_recursive = distilled_cuts
    min_capacity = capacities[1];
    new_max_capacity = total_wts - min_capacity 
    memento = Int[]
    no_soln = false
    capacities_recursive = [min_capacity, new_max_capacity]
    #println("Capacities starting ", capacities_recursive)
    for i in 1:recursion_levels
        (partition, cutsize, cutpoint) = two_way_linear_tree_sweep(tree_recursive, 
                                distilled_cuts_recursive,
                                hypergraph_recursive, 
                                capacities_recursive, 
                                2,
                                solns)
        #println("Recursive level ", i)
        #println("Capacities recursive  ", capacities_recursive)
        #println("Cutsize recorded ", cutsize)
        if cutsize >= 1e09
            no_soln = true
            break
        end
        blocks = zeros(Int, 2)
        for j in 1:length(partition) 
            blocks[partition[j]+1] += hypergraph_recursive.vwts[j]
        end
        (~, smaller_side) = findmin(blocks)
        recursive_total_vwts = 0
        for j in 1:hypergraph_recursive.num_vertices
            if partition[j] == smaller_side-1
                hypergraph_recursive.vwts[j] = 0
            else 
                recursive_total_vwts += hypergraph_recursive.vwts[j]
            end
        end
        distilled_cuts_recursive = distill_cuts_on_tree(hypergraph_recursive, 
                                                        fixed_vertices, 
                                                        T)
        capacities_recursive = [min_capacity, recursive_total_vwts-min_capacity]
        push!(memento, cutpoint)
        tree_recursive = deepcopy(T)
    end
    
    recursive_partition = -ones(Int, hgraph.num_vertices)
    cutsize = 1e09
    if no_soln == false
        for i in 1:length(memento)
            SimpleWeightedGraphs.rem_edge!(tree_recursive, memento[i], 
                                        distilled_cuts.pred[memento[i]])
        end
        comps = SimpleWeightedGraphs.connected_components(tree_recursive)
        recursive_partition = find_labels(comps, hgraph.num_vertices)
        (cutsize, balance) = golden_evaluator(hgraph, num_parts, recursive_partition)
        @debug "k-way tree sweep cut=$cutsize balance=$balance"
    else 
        @debug "k-way tree sweep found no valid balanced cut"
    end 
    return (recursive_partition, cutsize)
end

function generate_next_level(partition::Vector{Int},
                            hgraph::__hypergraph__,
                            T::SimpleWeightedGraphs.SimpleGraph,
                            original_capacities::Vector{Int},
                            capacities::Vector{Int},
                            num_parts::Int)
    balance = zeros(Int, num_parts)
    for i in eachindex(partition)
        balance[partition[i]+1] += hgraph.vwts[i]
    end
    (~, max_part) = findmax(balance)
    induced_hypergraph, clusters = extract_hypergraph(hgraph, partition, max_part-1)
    cluster_labels = findall(!iszero, clusters)
    new_min_capacity = capacities[1]
    new_max_capacity = sum(induced_hypergraph.vwts) - new_min_capacity
    new_capacities = [new_min_capacity, new_max_capacity]
    T_clustered = T[cluster_labels]
    comps = SimpleWeightedGraphs.connected_components(T_clustered)
    fixed = induced_hypergraph.fixed
    p1 = findall(x-> x == 0, fixed)
    p2 = findall(x-> x == 1, fixed)
    fixed_vertices = __pindex__(p1, p2)
    distilled_cuts_clustered = distill_cuts_on_tree(induced_hypergraph, 
                                                fixed_vertices, 
                                                T_clustered)
    return __recursive_parts__(induced_hypergraph,
                            T_clustered,
                            distilled_cuts_clustered,
                            new_capacities,
                            clusters)
end

# Process one (tree_type, eigenvector-subset) task: build the reweighted graph,
# construct the spanning tree, distill cuts, and produce candidate partitions
# from both the linear tree sweep and METIS-on-tree. Returns the candidate
# `[partition, cutsize]` entries for this task (1 or 2 of them).
# Number of randomized low-stretch (akpw) trees built per eigenvector subset.
# akpw is randomized, so distinct replicates yield distinct candidate trees and
# thus more diverse candidate partitions for the overlay/solve stage. MST trees
# (type 2) are deterministic, so only one replicate is built for them.
# Override with KSPECPART_LSST_TREES (set to 1 to recover the old behavior).
const NUM_LSST_TREES = parse(Int, get(ENV, "KSPECPART_LSST_TREES", "2"))

# Whether to also build one cut-aware-reweighted tree per (type, subset).
# Default off (opt-in via KSPECPART_CUTAWARE_TREES=1): did not improve results
# on the tested designs and adds tree/FM work. Code retained for experimentation.
const ENABLE_CUTAWARE_TREES = get(ENV, "KSPECPART_CUTAWARE_TREES", "0") == "1"

# Cap on the number of eigenvector subsets enumerated per tree type. The full
# enumeration is 2^d - 1 (d = embedding dimension); for the K-way path d = K-1,
# so without a cap large k explodes (e.g. k=64 -> 2^63 subsets). Below the cap
# we enumerate ALL subsets (preserving the original behavior for small d, e.g.
# 2-way eigvecs and k up to ~6); above it we use a bounded, informative set.
const MAX_TREE_SUBSETS = parse(Int, get(ENV, "KSPECPART_MAX_TREE_SUBSETS", "31"))

# Upper bound on the total number of tree candidates generated per refinement
# iteration. Each candidate is FM-refined (an external, per-candidate-expensive
# step at large k), so an unbounded count made k>=~6 runs impractical (~90+
# candidates). candidates ~ #subsets x (lsst_reps + 1).
const MAX_TREE_CANDIDATES = parse(Int, get(ENV, "KSPECPART_MAX_TREE_CANDIDATES", "24"))

# Auto-scale candidate generation with k. For small k (k <= 4) the natural
# subset count is tiny (2^(K-1)-1 <= 7) and the LSST replication is kept, so
# this is inert -- behavior is unchanged. For larger k (where the LDA embedding
# has ~K-1 dimensions and the subset enumerator would otherwise hit the full
# MAX_TREE_SUBSETS cap and be multiplied by the replicate count) it drops the
# LSST replicates to 1 and caps the subset count so the total candidate count
# stays near MAX_TREE_CANDIDATES.
function tree_generation_limits(num_parts::Int)
    reps = num_parts > 4 ? 1 : NUM_LSST_TREES
    max_subsets = max(1, min(MAX_TREE_SUBSETS, fld(MAX_TREE_CANDIDATES, reps + 1)))
    return max_subsets, reps
end

function tree_eigvec_subsets(dims::Vector{Int}, max_subsets::Int)
    d = length(dims)
    d <= 1 && return [copy(dims)]
    if (1 << min(d, 62)) - 1 <= max_subsets
        # All non-empty subsets (cheap for small d).
        subs = Vector{Vector{Int}}()
        for i in 1:d, s in Combinatorics.combinations(dims, i)
            push!(subs, s)
        end
        return subs
    end
    # Bounded set for large d: as many singletons as fit, plus the full set.
    subs = Vector{Vector{Int}}()
    for i in dims
        length(subs) >= max_subsets - 1 && break
        push!(subs, [i])
    end
    push!(subs, copy(dims))
    return subs
end

function process_tree_task(adj::SparseMatrixCSC,
                           X::Array{Float64},
                           hgraph::__hypergraph__,
                           fixed_vertices::__pindex__,
                           ub_factor::Real,
                           capacities::Vector{Int},
                           metis_path::String,
                           num_parts::Int,
                           seed::Int,
                           kway::Bool,
                           type::Int,
                           subset::Vector{Int},
                           rep::Int,
                           cutaware::Bool)
    results = Vector{Any}()
    # Seed this task's (task-local) RNG deterministically (including the
    # replicate index) so randomized akpw construction is both diverse across
    # replicates and reproducible regardless of how tasks are scheduled.
    Random.seed!(hash((seed, type, subset, rep, cutaware)))
    lst = type == 1
    X_thr = X[:, subset]
    if size(X_thr, 2) > 1 && type == 3
        return results
    end
    if type == 3
        X_thr = X_thr[:, 1]
    end
    @debug "tree task: eigenvectors=$(subset) type=$type cutaware=$cutaware"
    # Cut-aware reweighting blends the hyperedge-derived edge weights with the
    # embedding distance, biasing trees away from heavily-cut nets.
    clique_expansion = cutaware ? reweigh_graph_with_cuts(adj, hgraph, X_thr, lst) :
                                  reweigh_graph(adj, X_thr, lst)
    (tree, tree_matrix) = construct_tree(clique_expansion, X_thr, type)
    distilled_cuts = distill_cuts_on_tree(hgraph, fixed_vertices, tree)
    cutsize_tree = 1e09
    if kway == false
        (partition_tree, cutsize_tree, ~) = two_way_linear_tree_sweep(
            tree, distilled_cuts, hgraph, capacities, num_parts, 2)
    else
        (partition_tree, cutsize_tree) = k_way_linear_tree_sweep(
            tree, distilled_cuts, fixed_vertices, hgraph, capacities, num_parts, 2)
    end
    (partition_metis, cutsize_metis) = METIS_tree_partition(
        tree, distilled_cuts, hgraph, seed, metis_path, type, num_parts, ub_factor)
    push!(results, [partition_metis, cutsize_metis])
    if cutsize_tree < 1e09
        push!(results, [partition_tree, cutsize_tree])
    end
    return results
end

function tree_partition(adj::SparseMatrixCSC,
                        X::Array{Float64},
                        hgraph::__hypergraph__,
                        fixed_vertices::__pindex__,
                        ub_factor::Real,
                        capacities::Vector{Int},
                        metis_path::String,
                        num_parts::Int,
                        seed::Int,
                        kway::Bool)
    dims = Vector{Int}(1:size(X, 2))
    types = 2
    # The (tree_type, eigenvector-subset) tasks are independent; enumerate them
    # in the original order, then run them concurrently. Each task writes into
    # its own results slot (no shared push!), and METIS-on-tree now uses unique
    # temp filenames and cleans up after itself, so there is no shared state.
    max_subsets, lsst_reps = tree_generation_limits(num_parts)
    subsets = tree_eigvec_subsets(dims, max_subsets)
    tasks = Tuple{Int, Vector{Int}, Int, Bool}[]
    for type in 1:types
        # Randomized LSSTs (type 1) get multiple replicates; deterministic MSTs
        # (type 2) get a single one. Replicates are auto-scaled down at large k.
        nreps = type == 1 ? lsst_reps : 1
        for subset in subsets
            for rep in 1:nreps
                push!(tasks, (type, subset, rep, false))
            end
            # One extra cut-aware variant per (type, subset) for diversity.
            if ENABLE_CUTAWARE_TREES
                push!(tasks, (type, subset, 1, true))
            end
        end
    end
    task_results = Vector{Vector{Any}}(undef, length(tasks))
    @sync for k in 1:length(tasks)
        (type, subset, rep, cutaware) = tasks[k]
        Threads.@spawn begin
            task_results[k] = process_tree_task(adj, X, hgraph, fixed_vertices,
                ub_factor, capacities, metis_path, num_parts, seed, kway,
                type, subset, rep, cutaware)
        end
    end
    partitions = []
    # Direct spectral (k-means) clustering candidates, in addition to the
    # tree-sweep / METIS-on-tree candidates.
    append!(partitions, kmeans_candidates(X, hgraph, num_parts, capacities, seed))
    for res in task_results
        append!(partitions, res)
    end
    return partitions
end