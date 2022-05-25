@inline combine(x, y) = x

function GenerateEigenvectorAwareGraphAndTree(H::Hypergraph, X::Array{Float64}, lst::Bool)
    (~, nev) = size(X)
    offset = 0
    x = X[:, 1]
    i_tree = Vector{Int}()
    j_tree = Vector{Int}()

    v_ordering = sortperm(x)
    ordering_map = zeros(Int, H.n)

    for i in 1:length(v_ordering)
        ordering_map[v_ordering[i]] = i
    end
    
    for i in 1:H.e
        start_idx = H.eptr[i]
        end_idx = H.eptr[i+1]-1

        vertices = H.hedges[start_idx:end_idx]
        vertices_ordering = sort(ordering_map[vertices])

        for j in 1:length(vertices_ordering)-1
            v_start = v_ordering[vertices_ordering[j]]
            v_end = v_ordering[vertices_ordering[j+1]]
            push!(i_tree, v_start)
            push!(j_tree, v_end)
        end
    end

    A_tree = sparse(i_tree, j_tree, 1.0, H.n, H.n)
    A_tree = A_tree + A_tree'

    for ptr in 1:length(A_tree.colptr)-1
		row_len = A_tree.colptr[ptr+1]-A_tree.colptr[ptr]

		for row in 1:row_len
			row_ = A_tree.rowval[row+offset]
			path_distance = 0.0

			for d in 1:nev
				if lst == false
					path_distance += abs((X[row_, d] - X[ptr, d]))
				else
					if X[row_, d] - X[ptr, d] == 0.0
						path_distance += 1e09
					else
						path_distance += 1/abs(X[row_, d] - X[ptr, d])
					end
				end
            end
			A_tree.nzval[row+offset] = path_distance 
		end
		offset += row_len
	end

    r = akpw(A_tree)
	tree = LightGraphs.SimpleGraph(r)

    return tree
end