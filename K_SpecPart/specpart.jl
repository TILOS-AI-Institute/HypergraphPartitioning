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

function two_way_spectral_refine(hypergraph_file::String,
                                partition::Vector{Int},
                                hgraph::__hypergraph__,
                                metis_path::String,
                                ub_factor::Int;
                                num_parts::Int = 2,
                                eigen_vecs::Int = 1,
                                cycles::Int = 1,
                                refine_iters::Int = 2,
                                solver_iters::Int = 20,
                                best_solns::Int = 3,
                                seed::Int = 0
                                )
    line_log = repeat("=", 60)
    @info "$line_log"
    @info "**Solver parameters**"
    @info "$line_log"
    @info "Solver iterations $solver_iters"
    @info "Num vecs $eigen_vecs"
    @info "Tol 1e-40"
    processed_hg_name = hypergraph_file * ".processed"
    pre_refined_part = deepcopy(partition)
    pre_refined_cut, ~ = golden_evaluator(hgraph, num_parts, pre_refined_part)
    write_hypergraph(hgraph, processed_hg_name)
    specpart_refined_partition = partition
    adj_matrix = hypergraph2graph(hgraph, cycles)
    max_capacity = Int(ceil(sum(hgraph.vwts) * (50+ub_factor)/100))
    min_capacity = sum(hgraph.vwts) - max_capacity
    global_partitions = []
    global_cutsizes = Int[]
    for iter in 1:refine_iters
        @info "[specpart] iteration $iter"
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
        partitions = []
        cutsizes = Int[]
        cut_dictionary = Dict{Int, Int}()
        
        @sync Threads.@threads for i in 1:length(partitions_vec)
            partition_file = "tree_partition" * string(i) * ".part.2" 
            write_partition(partitions_vec[i][1], partition_file)
            triton_part_refine(triton_part_refiner_path, processed_hg_name, partition_file, num_parts, ub_factor, seed, i)
        end
        for i in 1:length(partitions_vec)
            partition_file = "tree_partition" * string(i) * ".part.2" 
            partitions_vec[i][1] = read_hint_file(partition_file)
            push!(partitions, partitions_vec[i][1])
        end
        rm_cmd = "rm -r *.part.2"
        run(`sh -c $rm_cmd`, wait=true)
        for i in 1:length(partitions)
            (cutsize, balance) = golden_evaluator(hgraph, num_parts, partitions[i])
            push!(cutsizes, cutsize)
            if haskey(cut_dictionary, cutsize) == false
                push!(cut_dictionary, cutsize => i)
            end
            @info "[specpart] Refined partition $i with cutsize $cutsize $balance" 
        end

        unique_partitions = []
        unique_cutsizes = Int[]
        unique_keys = collect(keys(cut_dictionary))
        for i in 1:length(unique_keys)
            key = unique_keys[i]
            partition_id = cut_dictionary[key]
            push!(unique_cutsizes, key)
            push!(unique_partitions, partitions[partition_id])
        end

        sorted_partition_ids = sortperm(unique_cutsizes)
        best_partitions = [] 
        solns_to_pick = min(best_solns, length(unique_cutsizes))
        for i in 1:solns_to_pick
            push!(best_partitions, unique_partitions[sorted_partition_ids[i]])
            @info "[specpart] partition picked with cutsize $(unique_cutsizes[sorted_partition_ids[i]])" 
        end
        push!(best_partitions, partition)
        hgraph_contracted, clusters = overlay(best_partitions, hgraph)
        @info "Running optimal attempt partitioning**"
        refined_partition = optimal_partitioner(hmetis_path, ilp_path, hgraph_contracted, num_parts, ub_factor)
        cutsize = golden_evaluator(hgraph_contracted, num_parts, refined_partition)
        @info "specpart cutsize recorded: $cutsize" 
        # projecting partition on clustered hypergraph to original hypergraph
        partition_projected = zeros(Int, hgraph.num_vertices)
        for i in 1:length(clusters)
            cid = clusters[i]
            partition_projected[i] = refined_partition[cid]
        end
        specpart_partition_name = processed_hg_name * ".specpart" * ".part.2"
        write_partition(partition_projected, specpart_partition_name)
        triton_part_refine(triton_part_refiner_path, processed_hg_name, specpart_partition_name, num_parts, ub_factor, seed, 0)
        partition_projected = read_hint_file(specpart_partition_name)
        cutsize = golden_evaluator(hgraph, num_parts, partition_projected)
        @info "specpart cutsize recorded: $cutsize" 
        specpart_refined_partition = partition_projected
        push!(global_partitions, partition_projected)
        push!(global_cutsizes, cutsize[1])
    end
    # global overlaying
    @info "[specpart] running final round of overlay"
    push!(global_partitions, partition)
    hgraph_contracted, clusters = overlay(global_partitions, hgraph)
    @info "Running optimal attempt partitioning**"
    specpart_partition = optimal_partitioner(hmetis_path, ilp_path, hgraph_contracted, num_parts, ub_factor)
    cutsize = golden_evaluator(hgraph_contracted, num_parts, specpart_partition)
    @info "[specpart] refined cutsize recorded: $cutsize" 
    partition_projected = zeros(Int, hgraph.num_vertices)
    for i in 1:length(clusters)
        cid = clusters[i]
        partition_projected[i] = specpart_partition[cid]
    end
    specpart_partition_name = processed_hg_name * ".specpart" * ".part.2"
    write_partition(partition_projected, specpart_partition_name)
    triton_part_refine(triton_part_refiner_path, processed_hg_name, specpart_partition_name, num_parts, ub_factor, seed, 0)
    partition_projected = read_hint_file(specpart_partition_name)
    cutsize = golden_evaluator(hgraph, num_parts, partition_projected)
    global_min_cut, idx = findmin(global_cutsizes)
    post_refined_cut = 0
    post_refined_part = Int[]
    if global_min_cut < cutsize[1]
        post_refined_cut = global_min_cut
        post_refined_part = global_partitions[idx]
    else
        post_refined_cut = cutsize[1]
        post_refined_part = partition_projected
    end 
    if pre_refined_cut < post_refined_cut
        @info "specpart cutsize recorded: $pre_refined_cut"
        return pre_refined_part, pre_refined_cut
    else 
        @info "specpart cutsize recorded: $post_refined_cut"   
        return post_refined_part, post_refined_cut
    end
end

function k_way_spectral_refine(hypergraph_file::String, 
                            partition::Vector{Int},
                            hgraph::__hypergraph__,
                            metis_path::String,
                            ub_factor::Int;
                            num_parts::Int = 2,
                            eigen_vecs::Int = 1,
                            cycles::Int = 1,
                            refine_iters::Int = 2,
                            solver_iters::Int = 20,
                            best_solns::Int = 3,
                            seed::Int = 0)
    line_log = repeat("=", 60)
    @info "$line_log"
    @info "**Solver parameters**"
    @info "$line_log"
    @info "Solver iterations $solver_iters"
    @info "Num vecs $eigen_vecs"
    @info "Tol 1e-40"
    specpart_refined_partition = partition
    processed_hg_name = hypergraph_file * ".processed"
    write_hypergraph(hgraph, processed_hg_name)
    adj_matrix = hypergraph2graph(hgraph, cycles)
    max_capacity = Int(ceil(sum(hgraph.vwts) * (50+ub_factor)/100))
    min_capacity = sum(hgraph.vwts) - max_capacity
    global_partitions = []
    global_cutsizes = Int[]
    for iter in 1:refine_iters
        @info "[specpart] iteration $iter"
        partition_list = [Vector{Int}() for i in 1:num_parts]
        for i in 1:length(partition)
            push!(partition_list[partition[i]+1], i)
        end
        hgraph_vec = []
        adj_mat_vec = []
        pindices_vec = []
        dims_ = Int(ceil(num_parts*(num_parts-1)*eigen_vecs)/2)
        for i in 1:num_parts
            side_0 = partition_list[i]
            side_1 = Int[]
            for j in 1:num_parts
                if i == j
                    continue
                end
                append!(side_1, partition_list[j])
                #side_0 = partition_list[i]
                #side_1 = partition_list[j]
                #fixed_vertices = __pindex__(side_0, side_1)
                #push!(pindices_vec, fixed_vertices)
                #push!(hgraph_vec, hgraph)
                #push!(adj_mat_vec, adj_matrix)
            end
            fixed_vertices = __pindex__(side_0, side_1)
            push!(pindices_vec, fixed_vertices)
            push!(hgraph_vec, hgraph)
            push!(adj_mat_vec, adj_matrix)
        end
        i_ptr = 0
        i_col_ptr = 0
        embedding_vec = []
        order_vec = []
        @sync Threads.@threads for i in 1:length(pindices_vec)
            embedding = solve_eigs(hgraph_vec[i], 
                                adj_mat_vec[i], 
                                pindices_vec[i], 
                                false, 
                                eigen_vecs, 
                                solver_iters,
                                epsilon=num_parts-1)
            push!(embedding_vec, embedding)
            push!(order_vec, i)
        end

        ordered_vec_idxs = sortperm(order_vec)
        concatenated_evec_matrix = embedding_vec[ordered_vec_idxs[1]]
        for i in 2:length(ordered_vec_idxs)
            embedding = embedding_vec[ordered_vec_idxs[i]]
            concatenated_evec_matrix = hcat(concatenated_evec_matrix, 
                                            embedding)
        end

        #evaluating clustering metrics
        #=pca_embedding = pca(concatenated_evec_matrix, partition)
        lda_embedding = lda(concatenated_evec_matrix, partition)
        reduced_evec_matrix = dimensionality_reduction(concatenated_evec_matrix, eigen_vecs, seed)
        reduced_evec_matrix = Array(reduced_evec_matrix')
        write_embedding_to_file(reduced_evec_matrix, "random_embedding.txt")
        write_embedding_to_file(pca_embedding, "pca_embedding.txt")
        write_embedding_to_file(lda_embedding, "lda_embedding.txt")
        write_labels_to_file(partition, "labels.txt")
        return concatenated_evec_matrix, partition=#

        reduced_evec_matrix = lda(concatenated_evec_matrix, partition)
        reduced_evec_matrix = Array(reduced_evec_matrix')
        # orthonormalization of the matrix
        #orthonormalized_evec_matrix = Matrix(qr(reduced_evec_matrix).Q)
        #filtered_evec_matrix = projection(concatenated_evec_matrix)
        #filtered_evec_matrix = concatenated_evec_matrix[:,1:eigen_vecs]
        max_capacity = Int(ceil(sum(hgraph.vwts) * ((100/num_parts) + ub_factor)/100))
        min_capacity = Int(ceil(sum(hgraph.vwts) * ((100/num_parts) - ub_factor)/100))
        fixed_vertices = __pindex__(Int[], Int[])
        partitions_vec = tree_partition(adj_mat_vec[1], 
                                        #concatenated_evec_matrix,
                                        reduced_evec_matrix,
                                        #filtered_evec_matrix, 
                                        hgraph, 
                                        fixed_vertices, 
                                        ub_factor,
                                        [min_capacity, max_capacity], 
                                        metis_path,
                                        num_parts,
                                        seed,
                                        true)
        #return partitions_vec
        partitions = []
        cutsizes = Int[]
        cut_dictionary = Dict{Int, Int}()

        #=Threads.@threads for i in 1:length(partitions_vec)
            partition_file = "tree_partition" * string(i) * ".part." * string(num_parts) 
            hypergraph_file_local = hypergraph_file * "." * string(i)
            hypergraph_file_copy = "cp $hypergraph_file $hypergraph_file_local"
            run(`sh -c $hypergraph_file_copy`, wait=true)
            write_partition(partitions_vec[i][1], partition_file)
            triton_part_refine(hypergraph_file_local, partition_file, num_parts, ub_factor, seed)
            partitions_vec[i][1] = read_hint_file(partition_file)
            push!(partitions, partitions_vec[i][1])
            rm_cmd = "rm -r $hypergraph_file_local"
            run(`sh -c $rm_cmd`, wait=true)
        end=#

        @sync Threads.@threads for i in 1:length(partitions_vec)
            partition_file = "tree_partition" * string(i) * ".part." * string(num_parts) 
            write_partition(partitions_vec[i][1], partition_file)
            triton_part_refine(triton_part_refiner_path, processed_hg_name, partition_file, num_parts, ub_factor, seed, i)
        end
        for i in 1:length(partitions_vec)
            partition_file = "tree_partition" * string(i) * ".part." * string(num_parts) 
            partitions_vec[i][1] = read_hint_file(partition_file)
            push!(partitions, partitions_vec[i][1])
        end

        rm_cmd = "rm -r *.part." * string(num_parts)
        run(`sh -c $rm_cmd`, wait=true)
        for i in 1:length(partitions)
            (cutsize, balance) = golden_evaluator(hgraph, num_parts, partitions[i])
            push!(cutsizes, cutsize)
            if haskey(cut_dictionary, cutsize) == false
                push!(cut_dictionary, cutsize => i)
            end
            @info "[specpart] Refined partition $i with cutsize $cutsize $balance" 
        end

        unique_partitions = []
        unique_cutsizes = Int[]
        unique_keys = collect(keys(cut_dictionary))
        for i in 1:length(unique_keys)
            key = unique_keys[i]
            partition_id = cut_dictionary[key]
            push!(unique_cutsizes, key)
            push!(unique_partitions, partitions[partition_id])
        end

        sorted_partition_ids = sortperm(unique_cutsizes)
        best_partitions = [] 
        solns_to_pick = min(best_solns, length(unique_cutsizes))
        for i in 1:solns_to_pick
            push!(best_partitions, unique_partitions[sorted_partition_ids[i]])
            @info "[specpart] partition picked with cutsize $(unique_cutsizes[sorted_partition_ids[i]])" 
        end
        push!(best_partitions, partition)
        hgraph_contracted, clusters = overlay(best_partitions, hgraph)
        #hgraph_contracted, clusters = iterative_overlay(best_partitions, hgraph)
        @info "Running optimal attempt partitioning**"
        refined_partition = optimal_partitioner(hmetis_path, ilp_path, hgraph_contracted, num_parts, ub_factor)
        cutsize = golden_evaluator(hgraph_contracted, num_parts, refined_partition)
        @info "specpart cutsize recorded: $cutsize" 
        # projecting partition on clustered hypergraph to original hypergraph
        partition_projected = zeros(Int, hgraph.num_vertices)
        for i in 1:length(clusters)
            cid = clusters[i]
            partition_projected[i] = refined_partition[cid]
        end
        #println("[debug], vertices ", hgraph.num_vertices, " , partition size ", length(partition_projected))
        specpart_partition_name = processed_hg_name * ".specpart" * ".part." * string(num_parts)
        write_partition(partition_projected, specpart_partition_name)
        triton_part_refine(triton_part_refiner_path, processed_hg_name, specpart_partition_name, num_parts, ub_factor, seed, 0)
        partition_projected = read_hint_file(specpart_partition_name)
        cutsize = golden_evaluator(hgraph, num_parts, partition_projected)
        @info "specpart cutsize recorded: $cutsize" 
        specpart_refined_partition = partition_projected
        push!(global_partitions, partition_projected)
        push!(global_cutsizes, cutsize[1])
    end
    # global overlaying
    @info "[specpart] running final round of overlay"
    push!(global_partitions, partition)
    hgraph_contracted, clusters = overlay(global_partitions, hgraph)
    @info "Running optimal attempt partitioning**"
    specpart_partition = optimal_partitioner(hmetis_path, ilp_path, hgraph_contracted, num_parts, ub_factor)
    cutsize = golden_evaluator(hgraph_contracted, num_parts, specpart_partition)
    @info "[specpart] refined cutsize recorded: $cutsize" 
    partition_projected = zeros(Int, hgraph.num_vertices)
    for i in 1:length(clusters)
        cid = clusters[i]
        partition_projected[i] = specpart_partition[cid]
    end
    pre_refined_part = partition_projected
    pre_refined_cut, ~ = golden_evaluator(hgraph, num_parts, pre_refined_part)
    specpart_partition_name = processed_hg_name * ".specpart" * ".part.2"
    write_partition(partition_projected, specpart_partition_name)
    triton_part_refine(triton_part_refiner_path, processed_hg_name, specpart_partition_name, num_parts, ub_factor, seed, 0)
    partition_projected = read_hint_file(specpart_partition_name)
    cutsize = golden_evaluator(hgraph, num_parts, partition_projected)
    global_min_cut, idx = findmin(global_cutsizes)
    final_cut = 0
    final_part = Int[]
    if cutsize[1] < global_min_cut
        final_cut = cutsize[1]
        final_part = partition_projected
    else
        final_cut = global_min_cut
        final_part = global_partitions[idx]
    end 
    cmd = "rm " * processed_hg_name
    run(`sh -c $cmd`, wait = true)
    post_refined_cut = final_cut
    post_refined_part = final_part
    if pre_refined_cut < post_refined_cut
        @info "specpart cutsize recorded: $pre_refined_cut" 
        return pre_refined_part, pre_refined_cut
    else    
        @info "specpart cutsize recorded: $post_refined_cut"
        return post_refined_part, post_refined_cut
    end
end

function specpart_run(hypergraph_file::String;
                    hypergraph_fixed_file::String = "",
                    hint_file::String = "",
                    imb::Int=2, 
                    num_parts::Int=2,
                    eigvecs::Int=2,
                    refine_iters::Int=2,
                    solver_iters::Int=40,
                    best_solns::Int=3,
                    ncycles = 1, 
                    seed::Int=0)
    BLAS.set_num_threads(Threads.nthreads())
    hg_split = split(hypergraph_file, "/")
    Random.seed!(seed)
    line_log = repeat("=", 60)
    @info "$line_log"
    @info "SpecPart K-way hypergraph partitioner**"
    @info "$line_log"
    @info "Hypergraph info**"
    @info "$line_log"
    hypergraph = read_hypergraph_file(hypergraph_file, hypergraph_fixed_file)
    @info "Num vertices $(hypergraph.num_vertices)"
    @info "Num hyperedges $(hypergraph.num_hyperedges)"
    @info "Fixed vertices $(length(findall(hypergraph.fixed .> -1)))" 
    (processed_hypergraph, 
    original_indices, 
    new_indices, 
    unused_indices) = isolate_islands(hypergraph)
    @info "$line_log"
    @info "Post processing**"
    @info "$line_log"
    @info "Num vertices $(processed_hypergraph.num_vertices)"
    @info "Num hyperedges $(processed_hypergraph.num_hyperedges)"
    partition = zeros(Int, hypergraph.num_vertices)
    if hint_file == ""
        @info "Running hmetis to generate a hint for SpecPart"
    else
        partition = read_hint_file(hint_file)
    end
    cutsize_pre_spec, ~ = golden_evaluator(hypergraph, num_parts, partition)
    @info "Hint partition cutsize $cutsize_pre_spec" 
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

