module SpecPart

using Shuffle 
using LinearAlgebra
using DataStructures
using SparseArrays
using Random
using Statistics
using LightGraphs
using SimpleWeightedGraphs
using Dates
using CRC32c

# Single source of truth for the version: the VERSION file at the repo root.
# Read once at load; fall back to "unknown" if the file is missing.
const KSPECPART_VERSION =
    try
        strip(read(joinpath(@__DIR__, "VERSION"), String))
    catch
        "unknown"
    end

include("logging_setup.jl")
include("definitions.jl")
include("dimensionality_reduction.jl")
include("run_triton_part_refiner.jl")
include("golden_evaluator.jl")
include("hypergraph.jl")
include("isolate_islands.jl")
include("io.jl")
include("graphification.jl")
include("embedding.jl")
include("tree_partition.jl")
include("overlay.jl")
include("projection.jl")

# ---------------------------------------------------------------------------
# Shared helpers for the two-way and K-way refinement loops. These factor out
# logic that was previously duplicated almost verbatim between the two.
# ---------------------------------------------------------------------------

# Tracks the best (lowest-cut) partition seen anywhere in the pipeline so the
# final result is the global best over every candidate, overlay result, and the
# hint. This makes the method monotone: adding more candidate solutions can only
# improve or tie the returned cut, never worsen it.
mutable struct __incumbent__
    cut::Int
    partition::Vector{Int}
end

function consider!(inc::__incumbent__, partition::Vector{Int}, cut::Int)
    if cut < inc.cut
        inc.cut = cut
        inc.partition = copy(partition)
    end
    return inc
end

# A partition is a valid k-way solution iff every label is in [0, num_parts),
# ALL k blocks are non-empty, and the per-block weight is within the balance
# tolerance (100/k + ub_factor)%. This rejects degenerate solutions -- most
# importantly the all-in-one-block partition, which has cut 0 and would
# otherwise be adopted as the global best.
function valid_kway_partition(partition::Vector{Int}, hgraph::__hypergraph__,
                              num_parts::Int, ub_factor::Real)
    blocks = zeros(Int, num_parts)
    @inbounds for i in 1:length(partition)
        b = partition[i]
        (b < 0 || b >= num_parts) && return false
        blocks[b+1] += hgraph.vwts[i]
    end
    all(>(0), blocks) || return false               # no empty block
    # Same capacity bound the rest of the pipeline targets (cf. check_balance).
    max_balance = Int(ceil((100/num_parts + ub_factor) * sum(hgraph.vwts) / 100))
    return all(b -> b <= max_balance, blocks)
end

# Offer a partition to the incumbent only if it is a valid k-way solution.
function consider_if_valid!(inc::__incumbent__, partition::Vector{Int}, cut::Int,
                            hgraph::__hypergraph__, num_parts::Int, ub_factor::Real)
    if valid_kway_partition(partition, hgraph, num_parts, ub_factor)
        consider!(inc, partition, cut)
    end
    return inc
end

# Best-effort repair of a candidate into a valid k-way partition: relabels any
# out-of-range labels, then moves vertices from the heaviest block to the
# lightest/empty blocks until every block is non-empty and the heaviest is
# within the balance tolerance. Idempotent for already-valid partitions (returns
# immediately). The moves are cut-blind -- cut quality is restored by the
# subsequent FM refinement; the goal here is only to turn a degenerate candidate
# (e.g. one with empty blocks) into a usable k-way partition instead of
# discarding it, and to give the external refiner a proper k-way input.
function repair_partition(partition::Vector{Int}, hgraph::__hypergraph__,
                          num_parts::Int, ub_factor::Real)
    n = length(partition)
    part = copy(partition)
    members = [Int[] for _ in 1:num_parts]
    w = zeros(Int, num_parts)
    @inbounds for v in 1:n
        b = part[v]
        if b < 0 || b >= num_parts
            b = 0
            part[v] = 0
        end
        push!(members[b+1], v)
        w[b+1] += hgraph.vwts[v]
    end
    max_balance = Int(ceil((100/num_parts + ub_factor) * sum(hgraph.vwts) / 100))
    for _ in 1:n
        heavy = argmax(w)
        light = argmin(w)
        (w[light] > 0 && w[heavy] <= max_balance) && break  # non-empty + balanced
        length(members[heavy]) <= 1 && break                # cannot move further
        v = pop!(members[heavy])
        part[v] = light - 1
        push!(members[light], v)
        w[heavy] -= hgraph.vwts[v]
        w[light] += hgraph.vwts[v]
    end
    return part
end

# FM-refine every candidate partition in parallel, then read the results back.
# Returns the refined partitions as a vector of label vectors.
function refine_candidates(partitions_vec, hgraph::__hypergraph__,
                           processed_hg_name::String,
                           num_parts::Int, ub_factor::Real, seed::Int)
    suffix = ".part." * string(num_parts)
    @sync Threads.@threads for i in 1:length(partitions_vec)
        partition_file = "tree_partition" * string(i) * suffix
        # Repair degenerate candidates (empty blocks / out-of-range labels) into
        # valid k-way partitions before refining, so every block is non-empty and
        # the external FM refiner receives a proper k-way input (it aborts at
        # large k when almost every block is empty).
        partitions_vec[i][1] = repair_partition(partitions_vec[i][1], hgraph,
                                                num_parts, ub_factor)
        write_partition(partitions_vec[i][1], partition_file)
        # If the external FM refiner fails anyway, keep the (repaired) unrefined
        # candidate rather than aborting the whole run.
        try
            triton_part_refine(triton_part_refiner_path, processed_hg_name,
                               partition_file, num_parts, ub_factor, seed, i)
        catch e
            @warn "FM refiner failed on candidate $i; keeping unrefined partition" exception=e
        end
    end
    partitions = Vector{Vector{Int}}()
    for i in 1:length(partitions_vec)
        partition_file = "tree_partition" * string(i) * suffix
        partitions_vec[i][1] = read_hint_file(partition_file)
        push!(partitions, partitions_vec[i][1])
    end
    run(`sh -c $("rm -r *" * suffix)`, wait=true)
    return partitions
end

# Evaluate candidates, offer each to the global incumbent, and return the
# `best_solns` lowest-cut (unique) partitions for the overlay. Empirically the
# overlay lift is dominated by feeding it LOW-CUT material, so selection is by
# cut (a diversity-aware variant was tried and consistently degraded the overlay
# because high-cut "diverse" candidates fragment the contracted hypergraph).
function select_best_partitions(partitions, hgraph::__hypergraph__,
                                num_parts::Int, best_solns::Int,
                                ub_factor::Real, inc::__incumbent__)
    cut_dictionary = Dict{Int, Int}()
    for i in 1:length(partitions)
        # Skip degenerate/invalid candidates so they are neither returned as the
        # global best nor fed (with an artificially low cut) into the overlay.
        valid_kway_partition(partitions[i], hgraph, num_parts, ub_factor) || continue
        (cutsize, balance) = golden_evaluator(hgraph, num_parts, partitions[i])
        consider!(inc, partitions[i], cutsize)
        if !haskey(cut_dictionary, cutsize)
            cut_dictionary[cutsize] = i
        end
        @debug "candidate $i refined cut=$cutsize balance=$balance"
    end
    unique_cutsizes = collect(keys(cut_dictionary))
    sorted_ids = sortperm(unique_cutsizes)
    best_partitions = Vector{Vector{Int}}()
    for i in 1:min(best_solns, length(unique_cutsizes))
        key = unique_cutsizes[sorted_ids[i]]
        push!(best_partitions, partitions[cut_dictionary[key]])
        @debug "selected candidate for overlay (cut=$key)"
    end
    return best_partitions
end

# Overlay-cluster `partitions`, solve the contracted hypergraph exactly,
# project back to the original vertices, FM-refine, and return the refined
# partition together with its (scalar) cut size.
function overlay_solve_project(partitions, hgraph::__hypergraph__,
                               processed_hg_name::String, num_parts::Int,
                               ub_factor::Real, seed::Int)
    hgraph_contracted, clusters = overlay(partitions, hgraph)
    @debug "overlay: $(hgraph_contracted.num_vertices) clusters; solving contracted hypergraph"
    refined_partition = optimal_partitioner(hmetis_path, ilp_path,
                                            hgraph_contracted, num_parts, ub_factor)
    @debug "contracted-solve cut=$(golden_evaluator(hgraph_contracted, num_parts, refined_partition))"
    partition_projected = zeros(Int, hgraph.num_vertices)
    for i in 1:length(clusters)
        partition_projected[i] = refined_partition[clusters[i]]
    end
    specpart_partition_name = processed_hg_name * ".specpart.part." * string(num_parts)
    write_partition(partition_projected, specpart_partition_name)
    try
        triton_part_refine(triton_part_refiner_path, processed_hg_name,
                           specpart_partition_name, num_parts, ub_factor, seed, 0)
    catch e
        @warn "FM refiner failed on overlay solution; keeping unrefined projection" exception=e
    end
    partition_projected = read_hint_file(specpart_partition_name)
    (cutsize, ~) = golden_evaluator(hgraph, num_parts, partition_projected)
    @debug "projected+refined cut=$cutsize"
    return partition_projected, cutsize
end

function two_way_spectral_refine(hypergraph_file::String,
                                partition::Vector{Int},
                                hgraph::__hypergraph__,
                                metis_path::String,
                                ub_factor::Real;
                                num_parts::Int = 2,
                                eigen_vecs::Int = 1,
                                cycles::Int = 1,
                                refine_iters::Int = 2,
                                solver_iters::Int = 20,
                                best_solns::Int = 3,
                                seed::Int = 0
                                )
    hgraph = remove_single_hyperedges(hgraph)
    processed_hg_name = hypergraph_file * ".processed"
    pre_refined_part = deepcopy(partition)
    pre_refined_cut, ~ = golden_evaluator(hgraph, num_parts, pre_refined_part)
    inc = __incumbent__(pre_refined_cut, copy(pre_refined_part))
    write_hypergraph(hgraph, processed_hg_name)
    specpart_refined_partition = partition
    adj_matrix = hypergraph2graph(hgraph, cycles)
    max_capacity = Int(ceil(sum(hgraph.vwts) * (50+ub_factor)/100))
    min_capacity = sum(hgraph.vwts) - max_capacity
    global_partitions = Vector{Vector{Int}}()
    global_cutsizes = Int[]
    for iter in 1:refine_iters
        @info "[iter $iter/$refine_iters] building and refining candidates"
        side_0 = findall(x-> x==0, specpart_refined_partition)
        side_1 = findall(x-> x==1, specpart_refined_partition)
        fixed_vertices = __pindex__(side_0, side_1)
        embedding = solve_eigs(hgraph, 
                            adj_matrix, 
                            fixed_vertices, 
                            false, 
                            eigen_vecs, 
                            solver_iters)
        fixed_vertices = __pindex__(Int[], Int[])
        partitions_vec = tree_partition(adj_matrix, 
                                        embedding, 
                                        hgraph, 
                                        fixed_vertices, 
                                        ub_factor,
                                        [min_capacity, max_capacity], 
                                        metis_path,
                                        num_parts,
                                        seed,
                                        false)
        partitions = refine_candidates(partitions_vec, hgraph, processed_hg_name,
                                       num_parts, ub_factor, seed)
        best_partitions = select_best_partitions(partitions, hgraph,
                                                 num_parts, best_solns, ub_factor, inc)
        push!(best_partitions, partition)
        partition_projected, cutsize = overlay_solve_project(best_partitions, hgraph,
                                            processed_hg_name, num_parts, ub_factor, seed)
        consider_if_valid!(inc, partition_projected, cutsize, hgraph, num_parts, ub_factor)
        @info "[iter $iter/$refine_iters] iteration cut = $cutsize | best so far = $(inc.cut)"
        specpart_refined_partition = partition_projected
        push!(global_partitions, partition_projected)
        push!(global_cutsizes, cutsize)
    end
    # global overlaying
    @info "[final] global overlay across all iterations"
    push!(global_partitions, partition)
    partition_projected, cutsize = overlay_solve_project(global_partitions, hgraph,
                                        processed_hg_name, num_parts, ub_factor, seed)
    consider_if_valid!(inc, partition_projected, cutsize, hgraph, num_parts, ub_factor)
    @debug "refine stage best cut=$(inc.cut)"
    return inc.partition, inc.cut
end

function k_way_spectral_refine(hypergraph_file::String, 
                            partition::Vector{Int},
                            hgraph::__hypergraph__,
                            metis_path::String,
                            ub_factor::Real;
                            num_parts::Int = 2,
                            eigen_vecs::Int = 1,
                            cycles::Int = 1,
                            refine_iters::Int = 2,
                            solver_iters::Int = 20,
                            best_solns::Int = 3,
                            seed::Int = 0)
    specpart_refined_partition = partition
    processed_hg_name = hypergraph_file * ".processed"
    hint_cut, ~ = golden_evaluator(hgraph, num_parts, partition)
    inc = __incumbent__(hint_cut, copy(partition))
    write_hypergraph(hgraph, processed_hg_name)
    adj_matrix = hypergraph2graph(hgraph, cycles)
    max_capacity = Int(ceil(sum(hgraph.vwts) * (50+ub_factor)/100))
    min_capacity = sum(hgraph.vwts) - max_capacity
    global_partitions = Vector{Vector{Int}}()
    global_cutsizes = Int[]
    for iter in 1:refine_iters
        @info "[iter $iter/$refine_iters] building and refining candidates"
        partition_list = [Vector{Int}() for i in 1:num_parts]
        for i in 1:length(partition)
            push!(partition_list[partition[i]+1], i)
        end
        # One eigenproblem per block: block i (side_0) vs. the rest (side_1).
        hgraph_vec = __hypergraph__[]
        adj_mat_vec = SparseMatrixCSC[]
        pindices_vec = __pindex__[]
        for i in 1:num_parts
            side_0 = partition_list[i]
            # Skip empty blocks: their supervision (bi_clique) term vanishes,
            # which yields a degenerate eigenproblem (common for large k with
            # sparse hints).
            isempty(side_0) && continue
            side_1 = Int[]
            for j in 1:num_parts
                j == i && continue
                append!(side_1, partition_list[j])
            end
            push!(pindices_vec, __pindex__(side_0, side_1))
            push!(hgraph_vec, hgraph)
            push!(adj_mat_vec, adj_matrix)
        end
        # Each block's eigenproblem writes into its own slot, so the threads
        # never mutate a shared container. The previous version used concurrent
        # `push!` into shared vectors (a data race) and then re-sorted them.
        n_problems = length(pindices_vec)
        embedding_vec = Vector{Matrix{Float64}}(undef, n_problems)
        # Parallelism here is at the block level, so drop BLAS to a single thread
        # inside the region to avoid n_problems x BLAS-threads oversubscription.
        blas_threads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        try
            @sync Threads.@threads for i in 1:n_problems
                embedding_vec[i] = solve_eigs(hgraph_vec[i],
                                    adj_mat_vec[i],
                                    pindices_vec[i],
                                    false,
                                    eigen_vecs,
                                    solver_iters,
                                    epsilon=num_parts-1)
            end
        finally
            BLAS.set_num_threads(blas_threads)
        end

        concatenated_evec_matrix = reduce(hcat, embedding_vec)

        # Supervised reduction of the per-block embeddings via multiclass LDA,
        # keeping up to (num_parts - 1) discriminative axes.
        reduced_evec_matrix = lda(concatenated_evec_matrix, partition, num_parts)
        reduced_evec_matrix = Array(reduced_evec_matrix')
        max_capacity = Int(ceil(sum(hgraph.vwts) * ((100/num_parts) + ub_factor)/100))
        min_capacity = Int(ceil(sum(hgraph.vwts) * ((100/num_parts) - ub_factor)/100))
        fixed_vertices = __pindex__(Int[], Int[])
        partitions_vec = tree_partition(adj_mat_vec[1],
                                        reduced_evec_matrix,
                                        hgraph,
                                        fixed_vertices,
                                        ub_factor,
                                        [min_capacity, max_capacity],
                                        metis_path,
                                        num_parts,
                                        seed,
                                        true)
        partitions = refine_candidates(partitions_vec, hgraph, processed_hg_name,
                                       num_parts, ub_factor, seed)
        best_partitions = select_best_partitions(partitions, hgraph,
                                                 num_parts, best_solns, ub_factor, inc)
        push!(best_partitions, partition)
        partition_projected, cutsize = overlay_solve_project(best_partitions, hgraph,
                                            processed_hg_name, num_parts, ub_factor, seed)
        consider_if_valid!(inc, partition_projected, cutsize, hgraph, num_parts, ub_factor)
        @info "[iter $iter/$refine_iters] iteration cut = $cutsize | best so far = $(inc.cut)"
        specpart_refined_partition = partition_projected
        push!(global_partitions, partition_projected)
        push!(global_cutsizes, cutsize)
    end
    # global overlaying
    @info "[final] global overlay across all iterations"
    push!(global_partitions, partition)
    partition_projected, cutsize = overlay_solve_project(global_partitions, hgraph,
                                        processed_hg_name, num_parts, ub_factor, seed)
    consider_if_valid!(inc, partition_projected, cutsize, hgraph, num_parts, ub_factor)
    cmd = "rm " * processed_hg_name
    run(`sh -c $cmd`, wait = true)
    @debug "refine stage best cut=$(inc.cut)"
    return inc.partition, inc.cut
end

function specpart_run(hypergraph_file::String;
                    hypergraph_fixed_file::String = "",
                    hint_file::String = "",
                    imb::Real=2, 
                    num_parts::Int=2,
                    eigvecs::Int=2,
                    refine_iters::Int=2,
                    solver_iters::Int=40,
                    best_solns::Int=5,
                    ncycles = 1, 
                    seed::Int=0)
    setup_logging()
    BLAS.set_num_threads(Threads.nthreads())
    Random.seed!(seed)
    run_start = time()
    @info LOG_RULE
    @info "K-SpecPart v$(KSPECPART_VERSION) : supervised spectral hypergraph partitioner"
    @info LOG_RULE
    hypergraph = read_hypergraph_file(hypergraph_file, hypergraph_fixed_file)
    log_kv("Input", basename(hypergraph_file))
    log_kv("Vertices", hypergraph.num_vertices)
    log_kv("Hyperedges", hypergraph.num_hyperedges)
    log_kv("Fixed verts", length(findall(hypergraph.fixed .> -1)))
    log_kv("Partitions", num_parts)
    log_kv("Imbalance", string(imb, "%"))
    # hMETIS/METIS treat imbalance RELATIVE to the ideal block weight (W/k). A
    # per-block slack `imb` larger than the ideal share (100/k) means a block may
    # exceed 2x the average, which lets the partitioners produce empty/degenerate
    # blocks at large k. Warn so a mis-specified imbalance is caught early.
    if num_parts > 1 && imb > 100 / num_parts
        rel = num_parts * imb / 100
        @warn @sprintf("imb=%.4g at k=%d is a very loose relative tolerance (heaviest block up to %.0f%% above average); this can yield empty/degenerate partitions. Consider a smaller imbalance.",
                       imb, num_parts, rel * 100)
    end
    log_kv("Parameters", "eigvecs=$eigvecs solver_iters=$solver_iters cycles=$ncycles " *
                         "best_solns=$best_solns refine_iters=$refine_iters")
    log_kv("Threads", Threads.nthreads())
    log_kv("Seed", seed)
    (processed_hypergraph,
    original_indices,
    new_indices,
    unused_indices) = isolate_islands(hypergraph)
    if processed_hypergraph.num_vertices != hypergraph.num_vertices
        log_kv("Largest CC", "$(processed_hypergraph.num_vertices) verts / " *
                             "$(processed_hypergraph.num_hyperedges) edges (islands removed)")
    end
    partition = zeros(Int, hypergraph.num_vertices)
    if hint_file == ""
        @warn "No hint file supplied; starting from a trivial partition"
    else
        partition = read_hint_file(hint_file)
    end
    cutsize_pre_spec, ~ = golden_evaluator(hypergraph, num_parts, partition)
    log_kv("Hint cut", cutsize_pre_spec)
    @info LOG_THIN
    processed_partition = -ones(Int, processed_hypergraph.num_vertices)
    process_hint(partition, new_indices, processed_partition) 
    refined_partition = similar(processed_partition)
    cutsize = 0
    if num_parts == 2
        refined_partition, cutsize = two_way_spectral_refine(hypergraph_file,
                                                    processed_partition, 
                                                    processed_hypergraph, 
                                                    metis_path,
                                                    imb, 
                                                    num_parts = num_parts,
                                                    eigen_vecs=eigvecs, 
                                                    refine_iters=refine_iters, 
                                                    solver_iters=solver_iters, 
                                                    cycles = ncycles,
                                                    best_solns=best_solns,
                                                    seed=seed)
    elseif num_parts > 2
        refined_partition, cutsize = k_way_spectral_refine(hypergraph_file,
                                                processed_partition, 
                                                processed_hypergraph, 
                                                metis_path,
                                                imb, 
                                                num_parts = num_parts,
                                                eigen_vecs=eigvecs, 
                                                refine_iters=refine_iters, 
                                                solver_iters=solver_iters, 
                                                cycles = ncycles,
                                                best_solns=best_solns,
                                                seed=seed)
    end

    elapsed = time() - run_start
    improvement = cutsize_pre_spec > 0 ?
        round(100 * (cutsize_pre_spec - cutsize) / cutsize_pre_spec, digits=2) : 0.0
    @info LOG_THIN
    log_kv("Hint cut", cutsize_pre_spec)
    log_kv("Final cut", cutsize)
    log_kv("Improvement", string(improvement, "%"))
    log_kv("Runtime", @sprintf("%.2f s", elapsed))
    @info LOG_RULE
    return refined_partition, cutsize

    output_ptn = zeros(Int, hypergraph.num_vertices)
    for i in 1:length(new_indices)
       output_ptn[i] = refined_partition[new_indices[i]]
    end
    for i in 1:length(unused_indices)
       output_ptn[i] = refined_partition[unused_indices[i]]
    end

    return output_ptn, cutsize
    #return refined_partition, cutsize
end
end

