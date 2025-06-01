@inline combine(x, y) = x
import Graphs

function ModifyExpanderWts(A::SparseMatrixCSC, eigvec::Array{Float64}, nev::Int, lst::Bool)
	n = size(A, 1)
	offset = 0
	A_ = deepcopy(A)
	g = SimpleWeightedGraph(A_)
	edge_weights = g.weights

	for ptr in 1:length(edge_weights.colptr)-1
		row_len = edge_weights.colptr[ptr+1]-edge_weights.colptr[ptr]

		for row in 1:row_len
			row_ = edge_weights.rowval[row+offset]
			distance = 0.0
			for d in 1:nev
				if lst == false
					distance += abs((eigvec[row_, d] - eigvec[ptr, d]))
				else
					if eigvec[row_, d] - eigvec[ptr, d] == 0.0
						distance += 1e09
					else
						distance += 1/abs(eigvec[row_, d] - eigvec[ptr, d])
					end
				end
			end
			edge_weights.nzval[row+offset] = distance
		end
		offset += row_len
	end
	
	g.weights = edge_weights

	return g
end

function GenEigenTree(g::SimpleWeightedGraph, A_wts::SparseMatrixCSC, eigvec::Array{Float64}, tree_type::Int)
	tree = SimpleWeightedGraphs.SimpleGraph()
	treeMatrix = spzeros(SimpleWeightedGraphs.nv(g), SimpleWeightedGraphs.nv(g))
	

	if tree_type == 1
		#r = akpw(g.weights)
		
		vtx_ids = Vector{Int}(1:SimpleWeightedGraphs.nv(g))
		eig_0 = eigvec[:, 1]
		vtx_ids = sortperm(eig_0)
		reverse_map = zeros(Int, length(vtx_ids))
		for i in 1:length(reverse_map)
			reverse_map[vtx_ids[i]] = i
		end
		i_g, j_g, w_g = findnz(g.weights)
		tree_reordered = sparse(vtx_ids[i_g], vtx_ids[j_g], w_g)
		tree_reordered_g = SimpleWeightedGraphs.SimpleWeightedGraph(tree_reordered)
		r = akpw(tree_reordered_g.weights)
		i_g, j_g, w_g = findnz(r)
		r = sparse(reverse_map[i_g], reverse_map[j_g], w_g)

		tree = SimpleWeightedGraphs.SimpleGraph(r)
		treeMatrix = r
		(i, j, v) = findnz(treeMatrix)
	elseif tree_type == 2
		#r = LightGraphs.kruskal_mst(g)
		#tree = SimpleGraph(nv(g))
		#r = SimpleWeightedGraphs.prim_mst(g)
		r = Graphs.prim_mst(g)
		i = Vector{Int}(undef, length(r))
		j = Vector{Int}(undef, length(r))
		v = Vector{Float64}(undef, length(r))

		@inbounds for idx in 1:length(r)
			vtxpair = r[idx]
			src_v = vtxpair.src
			dst_v = vtxpair.dst
			i[idx] = src_v
			j[idx] = dst_v
			v[idx] = g.weights[src_v, dst_v]
			#add_edge!(tree, r[idx].src, r[idx].dst)
		end

		tree = SimpleWeightedGraphs.SimpleGraph(r)
		treeMatrix = sparse(i, j, v, SimpleWeightedGraphs.nv(g), SimpleWeightedGraphs.nv(g))
		treeMatrix = treeMatrix + treeMatrix'
	elseif tree_type == 3
		r = modified_prim_mst(g, A_wts)
		i = Vector{Int}(undef, length(r))
		j = Vector{Int}(undef, length(r))
		v = Vector{Float64}(undef, length(r))

		@inbounds for idx in 1:length(r)
			vtxpair = r[idx]
			src_v = vtxpair.src
			dst_v = vtxpair.dst
			i[idx] = src_v
			j[idx] = dst_v
			v[idx] = g.weights[src_v, dst_v]
		end

		tree = SimpleWeightedGraphs.SimpleGraph(r)
		treeMatrix = sparse(i, j, v, nv(g), nv(g))
		treeMatrix = treeMatrix + treeMatrix'
	elseif tree_type == 4
		vertices = sortperm(eigvec)
		ii = vertices[1:end-1]
		jj = vertices[2:end]
		ww = abs.(eigvec[jj]-eigvec[ii])
		ww_z = findall(ww .== 0.0)
		ww[ww_z] .= 10^6
		treeMatrix = sparse(ii, jj, ww, SimpleWeightedGraphs.nv(g), SimpleWeightedGraphs.nv(g))
		treeMatrix = treeMatrix+treeMatrix'
		tree = SimpleWeightedGraphs.SimpleGraph(treeMatrix)
	else
		@warn "Please select correct Tree Type!! Exiting Partitioner!! "
	end

	return (tree, treeMatrix)
end

function GenerateShortestPathTrees(g::SimpleWeightedGraph, roots::Vector{Int})
	trees = Vector{SimpleGraph}()

	for i in 1:length(roots)
		ptr = 0
		ds = dijkstra_shortest_paths(g, roots[i])
		itree = zeros(Int, nv(g)-1)
		jtree = zeros(Int, nv(g)-1)

		for k in 1:length(ds.parents)

			if k == roots[i]
				continue
			end

			ptr += 1
			v = ds.parents[k]
			itree[ptr] = k
			jtree[ptr] = v
		end

		tree = sparse(itree, jtree, 1, nv(g), nv(g))
		tree = tree + tree'
		tree_graph = SimpleGraph(tree)
		push!(trees, tree_graph)
	end

	return trees
end

function GeneratePathTree(eigvec::Array{Float64})
	eigvec_ = eigvec[:, 1]
	n = length(eigvec_)
	nodes = sortperm(eigvec_)
	i = nodes[1:end-1]
	j = nodes[2:end]
	w = abs.(eigvec_[j]-eigvec_[i])
	w_z = findall(w .== 0.0)
	w[w_z] .= 10^6
	r = sparse(i, j, w, n, n)
	r = r+r'
    
	return (LightGraphs.SimpleGraph(r), r)
end

function GenEdgeWtMST(A::SparseMatrixCSC, lst::Bool)
	offset = 0
	n = size(A, 1)
	A_ = deepcopy(A)
	g = SimpleWeightedGraph(A_)
	edge_weights = g.weights

	for ptr in 1:length(edge_weights.colptr)-1
		row_len = edge_weights.colptr[ptr+1]-edge_weights.colptr[ptr]

		for row in 1:row_len
			val = edge_weights.nzval[row+offset]
			edge_weights.nzval[row+offset] = 1/val
		end
		offset += row_len
	end

	g.weights = edge_weights

	if lst==true
		r = akpw(g.weights)
		tree = LightGraphs.SimpleGraph(r)
		treeMatrix = r
		(i, j, v) = findnz(treeMatrix)
	else
		
		r = LightGraphs.prim_mst(g)
		i = Vector{Int}(undef, length(r))
		j = Vector{Int}(undef, length(r))
		v = Vector{Float64}(undef, length(r))

		@inbounds for idx in 1:length(r)
			vtxpair = r[idx]
			src_v = vtxpair.src
			dst_v = vtxpair.dst
			i[idx] = src_v
			j[idx] = dst_v
			v[idx] = A[src_v, dst_v]
		end

		tree = LightGraphs.SimpleGraph(r)
		treeMatrix = sparse(i, j, v, n, n)
		treeMatrix = treeMatrix + treeMatrix'
	end

	return (tree, treeMatrix)
end
