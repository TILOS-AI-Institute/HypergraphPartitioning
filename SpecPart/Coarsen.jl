function CoarsenHypergraph(hypergraph_c::Hypergraph_C, max_num_iter::Int, max_threshold::Int)
    record = Vector{Tuple{Int, Int}}()
    hypergraph = hypergraph_c.HG
    incidence_list = hypergraph_c.incidence_list
    hyperedges_pair_list = hypergraph_c.hyperedges_pair_list
    community = hypergraph_c.community
    fixed_part = hypergraph_c.fixed_part
    fixed_vertex_flag = hypergraph_c.fixed_vertex_flag

    memento = Vector{Vector{Tuple{Int, Int}}}()
    hedges = hypergraph.hedges
    hyperedges_weight = hypergraph.hwts
    vertices_weight = hypergraph.vwts
    vertices_flag = [false for i in 1:hypergraph.n]
    unvisited_vertices = Vector{Int}()
    num_remaining_vertices = 0
    num_fixed_vertices = 0
    hyperedges_flag = -ones(Int, hypergraph.e)
    hyperedges_parent = -ones(Int, hypergraph.e)

    if fixed_vertex_flag == true
        for i in 1:length(fixed_part)
            if fixed_part[i] > -1
                vertices_flag[i] = true
                num_fixed_vertices += 1
            else
                num_remaining_vertices += 1
            end
        end
    end

    for num_iter in 1:max_num_iter
        unvisited_vertices = findall(x-> x == false, vertices_flag)
        Shuffle.shuffle!(unvisited_vertices)

        if num_remaining_vertices == length(unvisited_vertices) && num_iter > 1
            break
        else
            num_remaining_vertices = length(unvisited_vertices)
        end
        
        num_remaining_hyperedges = 0
        
        for i in 1:length(hyperedges_flag)
            if hyperedges_flag[i] == -1
                num_remaining_hyperedges += 1
            end
        end

        log_syntax = repeat("*", 80)

        @info "$log_syntax"
        @info "[Coarsen] ITERATION : $num_iter with VERTICES : $num_remaining_vertices and HYPEREDGES : $num_remaining_hyperedges"

        record = Vector{Tuple{Int, Int}}()

        for i in 1:length(unvisited_vertices)
            u = unvisited_vertices[i]

            if vertices_flag[u] == true
                continue
            end

            hyperedges_u = incidence_list[u]
            cluster_score_map = Dict{Int, Float64}()

            for j in 1:length(hyperedges_u)
                hyperedge = hyperedges_u[j]

                if hyperedges_flag[hyperedge] > -1
                    continue
                end
                
                start_idx = hyperedges_pair_list[hyperedge, 1]
                end_idx = hyperedges_pair_list[hyperedge, 2]
                cluster_score = hyperedges_weight[hyperedge] * 1.0/(end_idx-start_idx)

                for idx in start_idx:end_idx
                    candidate = hedges[idx]

                    if u != candidate && vertices_flag[candidate] == false && 
                        community[u] == community[candidate] && 
                        vertices_weight[u] + vertices_weight[candidate] < max_threshold
                        
                        if haskey(cluster_score_map, candidate) == false
                            push!(cluster_score_map, candidate => cluster_score)
                        else
                            cluster_score_map[candidate] += cluster_score
                        end
                    end
                end
            end

            if isempty(cluster_score_map) == true
                continue
            end

            best_cluster_score = 0.0
            minimum_weight = -typemax(Float64)
            candidate = u

            for (v, score) in cluster_score_map
                if score > best_cluster_score
                    candidate = v
                    best_cluster_score = score
                    minimum_weight = vertices_weight[v]
                elseif score == best_cluster_score && minimum_weight < vertices_weight[v]
                    candidate = v
                    minimum_weight = vertices_weight[v]
                end
            end

            ContractVertex(u, candidate, hypergraph_c, hyperedges_flag, hyperedges_parent)
            vertices_flag[candidate] = true
            push!(record, (u, candidate))
        end

        push!(memento, record)
    end

    return (memento, vertices_flag, hyperedges_flag, hyperedges_parent)
end

function ContractVertex(u::Int, v::Int, hypergraph_c::Hypergraph_C, hyperedges_flag::Vector{Int}, hyperedges_parent::Vector{Int})
    hypergraph = hypergraph_c.HG
    incidence_list = hypergraph_c.incidence_list
    hyperedges_pair_list = hypergraph_c.hyperedges_pair_list
    hyperedges_hash = hypergraph_c.hyperedges_hash
    vertices_weight = hypergraph.vwts
    hyperedges_weight = hypergraph.hwts
    hedges = hypergraph.hedges
    hash_map = Dict{Int, Int}()
    vertices_weight[u] += vertices_weight[v]
    hyperedges_u = incidence_list[u]

    for i in 1:length(hyperedges_u)
        hyperedge = hyperedges_u[i]

        if hyperedges_flag[hyperedge] > -1
            continue
        end

        push!(hash_map, hyperedges_hash[hyperedge] => hyperedge)
    end

    hyperedges_v = incidence_list[v]

    for i in 1:length(hyperedges_v)
        hyperedge = hyperedges_v[i]

        if hyperedges_flag[hyperedge] > -1
            continue
        end

        start_idx = hyperedges_pair_list[hyperedge, 1]
        end_idx = hyperedges_pair_list[hyperedge, 2]
        u_idx = end_idx+1
        v_idx = end_idx+1

        for j in start_idx:end_idx
            if hedges[j] == u
                u_idx = j
            elseif hedges[j] == v
                v_idx = j
            end
        end

        if u_idx == end_idx+1
            # Vertex u is not in the hyperedge
            # Add vertex u to the bottom to make uncontraction easier
            hedges[v_idx] = hedges[end_idx]
            hedges[end_idx] = u 
            hyperedges_hash[hyperedge] += u*u - v*v
            push!(incidence_list[u], hyperedge)

            if haskey(hash_map, hyperedges_hash[hyperedge])
                hyperedges_flag[hyperedge] = v
                hyperedges_parent[hyperedge] = hash_map[hyperedges_hash[hyperedge]]
                hyperedges_weight[hash_map[hyperedges_hash[hyperedge]]] += hyperedges_weight[hyperedge]
            end
        else
            tmp = hedges[v_idx]
            hedges[v_idx] = hedges[end_idx]
            hedges[end_idx] = tmp
            hyperedges_pair_list[hyperedge, 2] -= 1
            hyperedges_hash[hyperedge] -= v*v
            
            if hyperedges_pair_list[hyperedge, 1] == hyperedges_pair_list[hyperedge, 2]
                hyperedges_flag[hyperedge] = v
                hyperedges_parent[hyperedge] = hyperedge
            elseif haskey(hash_map, hyperedges_hash[hyperedge])
                hyperedges_flag[hyperedge] = v
                hyperedges_parent[hyperedge] = hash_map[hyperedges_hash[hyperedge]]
                hyperedges_weight[hash_map[hyperedges_hash[hyperedge]]] += hyperedges_weight[hyperedge]
            end
        end
    end
end

#=
function CoarsenHypergraph(hypergraph_c::Hypergraph_C, max_num_iter::Int, max_threshold::Int)
    record = Vector{Tuple{Int, Int}}()
    hypergraph = hypergraph_c.HG
    incidence_list = hypergraph_c.incidence_list
    hyperedges_pair_list = hypergraph_c.hyperedges_pair_list
    community = hypergraph_c.community
    fixed_part = hypergraph_c.fixed_part
    fixed_vertex_flag = hypergraph_c.fixed_vertex_flag

    memento = Vector{Vector{Tuple{Int, Int}}}()
    hedges = hypergraph.hedges
    hyperedges_weight = hypergraph.hwts
    vertices_weight = hypergraph.vwts
    vertices_flag = [false for i in 1:hypergraph.n]
    unvisited_vertices = Vector{Int}()
    num_remaining_vertices = 0
    num_fixed_vertices = 0
    hyperedges_flag = -ones(Int, hypergraph.e)

    if fixed_vertex_flag == true
        for i in 1:length(fixed_part)
            if fixed_part[i] > -1
                vertices_flag[i] = true
                num_fixed_vertices += 1
            else
                num_remaining_vertices += 1
            end
        end
    end

    for num_iter in 1:max_num_iter
        unvisited_vertices = findall(x-> x == false, vertices_flag)
        Shuffle.shuffle!(unvisited_vertices)

        if num_remaining_vertices == length(unvisited_vertices)
            break
        else
            num_remaining_vertices = length(unvisited_vertices)
        end
        
        num_remaining_hyperedges = 0
        
        for i in 1:length(hyperedges_flag)
            if hyperedges_flag[i] == -1
                num_remaining_hyperedges += 1
            end
        end

        log_syntax = repeat("*", 80)

        @info "$log_syntax"
        @info "[Coarsen] ITERATION : $num_iter ;; VERTICES : $num_remaining_vertices :: HYPEREDGES : $num_remaining_hyperedges"

        record = Vector{Tuple{Int, Int}}()

        for i in 1:length(unvisited_vertices)
            u = unvisited_vertices[i]

            if vertices_flag[u] == true
                continue
            end

            hyperedges_u = incidence_list[u]
            cluster_score_map = Dict{Int, Float64}()

            for j in 1:length(hyperedges_u)
                hyperedge = hyperedges_u[j]

                if hyperedges_flag[hyperedge] > -1
                    continue
                end
                
                start_idx = hyperedges_pair_list[hyperedge, 1]
                end_idx = hyperedges_pair_list[hyperedge, 2]
                cluster_score = hyperedges_weight[hyperedge] * 1.0/(end_idx-start_idx)

                for idx in start_idx:end_idx
                    candidate = hedges[idx]

                    if u != candidate && vertices_flag[candidate] == false && 
                        community[u] == community[candidate] && 
                        vertices_weight[u] + vertices_weight[candidate] < max_threshold
                        
                        if haskey(cluster_score_map, candidate) == false
                            push!(cluster_score_map, candidate => cluster_score)
                        else
                            cluster_score_map[candidate] += cluster_score
                        end
                    end
                end
            end

            if isempty(cluster_score_map) == true
                continue
            end

            best_cluster_score = 0.0
            minimum_weight = -typemax(Float64)
            candidate = u

            for (v, score) in cluster_score_map
                if score > best_cluster_score
                    candidate = v
                    best_cluster_score = score
                    minimum_weight = vertices_weight[v]
                elseif score == best_cluster_score && minimum_weight < vertices_weight[v]
                    candidate = v
                    minimum_weight = vertices_weight[v]
                end
            end

            ContractVertex(u, candidate, hypergraph_c, hyperedges_flag)
            vertices_flag[candidate] = true
            push!(record, (u, candidate))
        end

        push!(memento, record)
    end

    return (memento, vertices_flag, hyperedges_flag)
end

function ContractVertex(u::Int, v::Int, hypergraph_c::Hypergraph_C, hyperedges_flag::Vector{Int})
    hypergraph = hypergraph_c.HG
    incidence_list = hypergraph_c.incidence_list
    hyperedges_pair_list = hypergraph_c.hyperedges_pair_list
    hyperedges_hash = hypergraph_c.hyperedges_hash
    vertices_weight = hypergraph.vwts
    hyperedges_weight = hypergraph.hwts
    hedges = hypergraph.hedges

    vertices_weight[u] += vertices_weight[v]
    hyperedges_v = incidence_list[v]

    for i in 1:length(hyperedges_v)
        hyperedge = hyperedges_v[i]

        if hyperedges_flag[hyperedge] > -1
            continue
        end

        start_idx = hyperedges_pair_list[hyperedge, 1]
        end_idx = hyperedges_pair_list[hyperedge, 2]

        u_idx = end_idx+1
        v_idx = end_idx+1

        for j in start_idx:end_idx
            if hedges[j] == u
                u_idx = j
            elseif hedges[j] == v
                v_idx = j
            end
        end

        if u_idx == end_idx+1
            # Vertex u is not in the hyperedge
            # Add vertex u to the bottom to make uncontraction easier
            hedges[v_idx] = hedges[end_idx]
            hedges[end_idx] = u 
            hyperedges_hash[hyperedge] += u*u - v*v
            push!(incidence_list[u], hyperedge)
        else
            tmp = deepcopy(hedges[v_idx])
            hedges[v_idx] = hedges[end_idx]
            hedges[end_idx] = tmp
            hyperedges_pair_list[hyperedge, 2]-= 1
            hyperedges_hash[hyperedge] -= v*v
            # Check if single hyperdge
            if hyperedges_pair_list[hyperedge, 1]  == hyperedges_pair_list[hyperedge, 2]
                hyperedges_flag[hyperedge] = hyperedge
            end
        end
    end

    # Detect parallel hyperedges
    hash_map = Dict{Int, Int}()
    net_map = Dict{Int, Int}()

    for i in 1:length(incidence_list[u])
        hyperedge = incidence_list[u][i]

        if hyperedges_flag[hyperedge] > -1 
            continue
        end

        if haskey(hash_map, hyperedges_hash[hyperedge]) == false
            push!(hash_map, hyperedges_hash[hyperedge] => hyperedge)
        else
            push!(net_map, hyperedge => hash_map[hyperedges_hash[hyperedge]])
        end
    end

    for (hyperedge, match_hyperedge) in net_map
        hyperedges_weight[match_hyperedge] += hyperedges_weight[hyperedge]
        hyperedges_flag[hyperedge] = match_hyperedge
    end
end=#