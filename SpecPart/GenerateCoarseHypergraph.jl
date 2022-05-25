function GenerateCoarseHypergraph(h_c::Hypergraph_C, vertices_flag::Vector{Bool}, hyperedges_flag::Vector{Int})
    vertex_id = 0
    hypergraph = h_c.HG
    fixed_vertex_flag = h_c.fixed_vertex_flag
    fixed_part = h_c.fixed_part
    vertex_weights = hypergraph.vwts
    hyperedge_weights = hypergraph.hwts
    vertex_weights_cc = zeros(Int, hypergraph.n)
    vertex_map = zeros(Int, hypergraph.n)

    for i in 1:length(vertices_flag)
        if (fixed_vertex_flag == true && fixed_part[i] > -1) || vertices_flag[i] == false
            vertex_id += 1
            vertex_weights_cc[vertex_id] = vertex_weights[i]
            vertex_map[i] = vertex_id
        end
    end
    
    hyperedge_id = 0
    hyperedge_weights_cc = similar(hyperedge_weights)
    hedges = hypergraph.hedges
    hedges_cc = similar(hedges)
    eptr = hypergraph.eptr
    eptr_cc = similar(eptr)
    offset = 0
    cc_k = 0

    for i in 1:length(hyperedges_flag)
        if hyperedges_flag[i] == -1
            hyperedge_id += 1
            hyperedge_weights_cc[hyperedge_id] = hyperedge_weights[i]
            start_idx = eptr[i]
            end_idx = eptr[i+1]-1
            hsize = 0

            for j in start_idx:end_idx
                v = hedges[j]

                if vertices_flag[v] == false
                    cc_k += 1
                    hedges_cc[cc_k] = vertex_map[v]
                    hsize += 1
                end
            end

            eptr_cc[hyperedge_id] = offset+1
            offset += hsize
        end
    end

    eptr_cc[hyperedge_id+1] = offset+1
    fixed_part_cc = zeros(Int, vertex_id)

    for i in 1:length(vertices_flag)
        if vertices_flag[i] == false
            vertex_mapped = vertex_map[i]
            fixed_part_cc[vertex_mapped] = fixed_part[i]
        end
    end

    return Hypergraph(vertex_id, hyperedge_id, hedges_cc[1:cc_k], eptr_cc[1:hyperedge_id+1], vertex_weights_cc[1:vertex_id], hyperedge_weights_cc[1:hyperedge_id]), fixed_part_cc 
end
