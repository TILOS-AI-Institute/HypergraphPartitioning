function UncoarsenHypergraph(h_c::Hypergraph_C, memento::Vector{Vector{Tuple{Int, Int}}}, vertices_flag::Vector{Bool}, hyperedges_flag::Vector{Int}, hyperedges_parent::Vector{Int}, part_vector::Vector{Int}, part_area::Vector{Int}, max_capacity::Int, min_capacity::Int)
    log_syntax = repeat("*", 80)
    total_vertices = length(findall(vertices_flag .== false))
    total_hyperedges = length(findall(hyperedges_flag .== -1))
    iter = 1

    hypergraph_c = h_c
    
    @info "$log_syntax"
    @info "[Uncoarsen] ITERATION : $iter with VERTICES : $total_vertices and HYPEREDGES : $total_hyperedges"

    for i in length(memento)-1:-1:1
        iter += 1
        cluster_records = memento[i]

        for j in length(cluster_records):-1:1
            (total_vertices, total_hyperedges) = UncontractVertex(cluster_records[j][1], cluster_records[j][2], vertices_flag, hyperedges_flag, hyperedges_parent, h_c, total_vertices, total_hyperedges)
        end

        @info "$log_syntax"
        @info "[Uncoarsen] ITERATION : $iter with VERTICES : $total_vertices and HYPEREDGES : $total_hyperedges"

        (hypergraph_clustered, fixed_part_clustered) = GenerateCoarseHypergraph(hypergraph_c, vertices_flag, hyperedges_flag)
        incidence_struct_clustered = HypergraphToIncidence(hypergraph_clustered)
        incidence_list_clustered = GenerateIncidenceList(incidence_struct_clustered)
        hyperedge_pair_list_clustered = GenerateHypergraphPairList(hypergraph_clustered)
        hyperedges_hash_clustered = GenerateHyperedgesHash(hypergraph_clustered)
        hypergraph_clustered_c = Hypergraph_C(hypergraph_clustered, incidence_list_clustered, hyperedge_pair_list_clustered, 
                    ones(Int, hypergraph_clustered.n), fixed_part_clustered, false, hyperedges_hash_clustered)

        (part_vector, global_min_cut_size, global_part_area) = FM(hypergraph_clustered_c, incidence_struct_clustered, part_vector, part_area, 
                                                                    max_capacity, min_capacity)

        @info "[CUT] FM CUT IS $global_min_cut_size :: $(global_part_area[1]) || $(global_part_area[2])"
        break
    end
end

function UncontractVertex(u::Int, v::Int, vertices_flag::Vector{Bool}, hyperedges_flag::Vector{Int}, hyperedges_parent::Vector{Int}, hypergraph_c::Hypergraph_C, total_vertices::Int, total_hyperedges::Int)
    hypergraph = hypergraph_c.HG
    incidence_list = hypergraph_c.incidence_list
    hyperedges_pair_list = hypergraph_c.hyperedges_pair_list
    hyperedges_hash = hypergraph_c.hyperedges_hash
    vertices_weight = hypergraph.vwts
    hyperedges_weight = hypergraph.hwts
    vertices_weight[u] -= vertices_weight[v]
    vertices_flag[v] = false
    total_vertices += 1
    hyperedges_v = incidence_list[v]
    k = length(hypergraph.hedges)

    for i in length(hyperedges_v):-1:1
        hyperedge = hyperedges_v[i]
        flag = hyperedges_flag[hyperedge] == -1 || hyperedges_flag[hyperedge] == v

        if flag == true
            if hyperedges_flag[hyperedge] == v && hyperedges_parent[hyperedge] != hyperedge
                hyperedges_weight[hyperedges_parent[hyperedge]] -= hyperedges_weight[hyperedge]
                total_hyperedges += 1
            elseif hyperedges_flag[hyperedge] == v && hyperedges_parent[hyperedge] == hyperedge
                total_hyperedges += 1
            end

            hyperedges_flag[hyperedge] = -1
            end_idx = hyperedges_pair_list[hyperedge, 2]

            if end_idx < k
                if (hypergraph.hedges[end_idx+1]) == v && ((end_idx+1) < hypergraph.eptr[hyperedge+1])
                    hyperedges_pair_list[hyperedge, 2] += 1
                    hyperedges_hash[hyperedge] += v*v
                else
                    h_start_idx = hyperedges_pair_list[hyperedge, 1]
                    h_end_idx  = hyperedges_pair_list[hyperedge, 2]

                    for j in h_end_idx:-1:h_start_idx
                        if hypergraph.hedges[j] == u
                            hypergraph.hedges[j] = v
                            break
                        end
                    end

                    hyperedges_hash[hyperedge] += v*v - u*u 
                    popat!(incidence_list[u], length(incidence_list[u]))
                end
            end
        end
    end

    return (total_vertices, total_hyperedges)
end


#=
function UncontractVertex(u::Int, v::Int, vertices_flag::Vector{Bool}, hyperedges_flag::Vector{Int}, hyperedges_parent::Vector{Int}, hypergraph_c::Hypergraph_C, total_vertices::Int, total_hyperedges::Int)
    hypergraph = hypergraph_c.HG
    incidence_list = hypergraph_c.incidence_list
    hyperedges_pair_list = hypergraph_c.hyperedges_pair_list
    hyperedges_hash = hypergraph_c.hyperedges_hash
    vertices_weight = hypergraph.vwts
    hyperedges_weight = hypergraph.hwts
    vertices_weight[u] -= vertices_weight[v]
    hyperedges_v = incidence_list[v]
    vertices_flag[v] = false
    total_vertices += 1
    
    #@info "[UC Debug] Hyperedges of $u : $(incidence_list[u])"
    #@info "[UC Debug] Hyperedges of $v : $(incidence_list[v])"
    #@info "[UC Debug] Common Hyperedges : $(intersect(incidence_list[u],incidence_list[v]))"

    for i in length(hyperedges_v):-1:1
        hyperedge = hyperedges_v[i]
        end_idx = hyperedges_pair_list[hyperedge, 2]

        #@info "[UC Debug] UNCONTRACTING: $u :: $v  WITH FLAGS $(vertices_flag[u]) :: $(vertices_flag[v]) "
        #@info "[UC Debug] Hyperedge id $hyperedge with flag $(hyperedges_flag[hyperedge])"
        #@info "[UC Debug] Original hyperedge: $(hypergraph.hedges[hypergraph.eptr[hyperedge] : hypergraph.eptr[hyperedge+1]-1])"
        #@info "[UC Debug] Before Contracting $(hypergraph.hedges[hyperedges_pair_list[hyperedge, 1]:hyperedges_pair_list[hyperedge, 2]])"
        #@info "[UC Debug] Check conditions $(hypergraph.hedges[end_idx+1]) :: $v :: $(end_idx+1) :: $(hypergraph.eptr[hyperedge+1])"

        if (hypergraph.hedges[end_idx+1]) == v && ((end_idx+1) < hypergraph.eptr[hyperedge+1])
            hyperedges_pair_list[hyperedge, 2] += 1
            hyperedges_hash[hyperedge] += v*v

            if hyperedges_flag[hyperedge] == v 
                hyperedges_flag[hyperedge] = -1
                total_hyperedges += 1
            end
        else
            if hyperedges_flag[hyperedge] == v || hyperedges_flag[hyperedge] == -1
                h_start_idx = hyperedges_pair_list[hyperedge, 1]
                h_end_idx  = hyperedges_pair_list[hyperedge, 2]
                u_flag = false

                for j in h_end_idx:-1:h_start_idx
                    if hypergraph.hedges[j] == u
                        hypergraph.hedges[j] = v
                        u_flag = true
                        break
                    end
                end

                if u_flag == true
                    hyperedges_hash[hyperedge] += v*v - u*u 

                    if hyperedges_flag[hyperedge] > -1
                        hyperedges_weight[hyperedges_flag[hyperedge]] -= hyperedges_weight[hyperedge]
                        hyperedges_flag[hyperedge] = -1
                        total_hyperedges += 1
                        x = popat!(incidence_list[u], length(incidence_list[u]))
                        #@info "[UC Debug] Popping $x"
                    end
                end
            else
                #=hyperedges_weight[hyperedges_flag[hyperedge]] -= hyperedges_weight[hyperedge]
                hyperedges_flag[hyperedge] = -1
                total_hyperedges += 1=#
            end
        end

        #@info "[UC Debug] After Contracting $(hypergraph.hedges[hyperedges_pair_list[hyperedge, 1]:hyperedges_pair_list[hyperedge, 2]]) with flag $(hyperedges_flag[hyperedge])"
        #@info "***********************************************"

    end

    #@info "[UC Debug] Hyperedges of $u : $(incidence_list[u])"
    #@info "[UC Debug] Hyperedges of $v : $(incidence_list[v])"

    return (total_vertices, total_hyperedges)
end
=#