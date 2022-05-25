function GenerateTwoPinHyperedges(hypergraph_c::Hypergraph_C, max_iters::Int)
    hypergraph = hypergraph_c.HG
    hyperedges_pair_list - hypergraph_c.hyperedges_pair_list
    incidence_list = hypergraph_c.incidence_list
    hsizes = hypergraph.eptr[2:end] - hypergraph.eptr[1:end-1]
    multi_pin_hyperedges = findall(x->x > 2, hsizes)
    total_two_pin_hyperedges = hypergraph.e - multi_pin_hyperedges
    hyperedges_ordering = sortperm(multi_pin_hyperedges)
    clusters = Vector{Int}(1:hypergraph.n)
    hyperedges_flag = -ones(Int, hypergraph.e)
    cluster_hyperedge_list = [Vector{Int}() for i in 1:hypergraph.e]
    cluster_id = 1

    for iter in 1:max_iters
        recheck_hyperedges = Vector{Int}()

        for i in 1:length(hyperedges_ordering)
            hyperedge = hyperedges_ordering[i]

            if hyperedges_flag[hyperedge] > -1
                continue
            end
            
            cant_cluster_all_flag = false
            
            start_idx = hyperedges_pair_list[hyperedge, 1]
            end_idx = hyperedges_pair_list[hyperedge, 2]

            for j in start_idx:end_idx
                vtx = hypergraph.hedges[j]

                if vertices_flag[vtx] > 0
                    cant_cluster_all_flag = true
                    push!(cluster_hyperedge_list[hyperedge], vtx)
                    break
                end

                vertices_flag[vtx] = 
                clusters[vtx] = cluster_id
            end

            if cant_cluster_all_flag == true
                push!(recheck_hyperedges, hyperedge)
            else
                hyperedges_flag[hyperedge] = cluster_id
                cluster_id += 1
            end
        end

        for i in 1:length(recheck_hyperedges)
            hyperedge = recheck_hyperedges[i]
            start_idx = hyperedges_pair_list[hyperedge, 1]
            end_idx = hyperedges_pair_list[hyperedge, 2]

            for j in start_idx:end_idx
                vtx = hypergraph.hedges[j]

                if vertices_flag[vtx] > 0 
                    continue
                end

                clusters[vtx] = cluster_id
            end

            cluster_id += 1
        end
    end
end