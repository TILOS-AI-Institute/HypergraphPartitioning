function findLabels(clusters::AbstractArray, n::Int)
	label = zeros(Int, n)

	for i in 1:length(clusters)
		for j in 1:length(clusters[i])
			label[clusters[i][j]] = i
		end
	end
    
	return label
end

function incidentEdges(H::Hypergraph, B::Incidence, S::Vector{Int}, hyperedges_flag::Vector{Int})
    hedges = B.hedges
	eptr = B.eptr
	e = H.e
	T = 0
	
    for j in 1:length(S)
        k = S[j]
        T += eptr[k+1]-eptr[k]
    end    

    E = zeros(Int, T)
	r = 0
	
    for j in 1:length(S)
        k = S[j]
        hedges_j = hedges[eptr[k]:(eptr[k+1]-1)]

        for i in 1:length(hedges_j)
            if hyperedges_flag[i] > 0
                r += 1
                E[r] = hedges_j[i]
            end
        end 
    end      
    
	unique!(E)
	return E
end

function incidentNodes(H::Hypergraph, E::Vector{Int})
    n = H.n
    eptr = H.eptr
    hedges = H.hedges
    
	T = 0
	
    for j in 1:length(E)
        k = E[j]
        T += eptr[k+1]-eptr[k]
    end

    S = zeros(Int, T)
	r = 1
	
    for j in 1:length(E)
        k = E[j]
        L = eptr[k+1]-eptr[k]
        S[r:r+L-1] = hedges[eptr[k]:(eptr[k+1]-1)]
        r += L 
    end      
    
    unique!(S)
	return S
end

function computeEdgeCuts(H::Hypergraph, lcaInstance::LCA, E::Vector{Int}, hyperedges_flag::Vector{Int})
	w_ = H.hwts
	we = 0
	lca_R = lcaInstance.fts
	lca_ET = lcaInstance.eulerTour
	lca_LV = lcaInstance.levelVec
	lca_M = lcaInstance.rmqSparseTable
	ifts = lcaInstance.ifts
	fts = lcaInstance.fts
	n = H.n
	M = H.e
	v_ = H.hedges
	loc_ = H.eptr
	edgeDiff = zeros(Int, n)
	edgeTerminators = zeros(Int, M)
	heap = zeros(Int, max(maximum(loc_[2:end]-loc_[1:end-1]), 1))
	heapIndex = 0

	@inbounds for jj in 1:length(E) 
		j = E[jj]

        if hyperedges_flag[j] == 0
            continue
        end

		we = w_[j]
		l_edge = loc_[j+1]-loc_[j]
		
		if l_edge == 2   
			u = v_[loc_[j]]
			v = v_[loc_[j]+1]
			w = lcaQuery(u, v, lca_LV, lca_ET, lca_R, lca_M)
			edgeTerminators[j] = w 

			if w == u
				edgeDiff[u] -= we
				edgeDiff[v] += we
			elseif w == v
				edgeDiff[v] -= we
				edgeDiff[u] += we
			else
				edgeDiff[v] += we
				edgeDiff[u] += we
				edgeDiff[w] -= 2*we
			end
		else
			for k in 0:l_edge-1
				u = v_[loc_[j]+k]
				edgeDiff[u] += we
				u = fts[u]
				x = u
				heapIndex += 1
				j0 = heapIndex

				while j0 > 1
					pj = Int(floor(j0/2))

					if heap[pj] > x
						heap[j0] = heap[pj]
					else
						break     
					end

					j0 = pj
				end

				heap[j0] = x
			end 

			while heapIndex > 1
				u = heap[1]
				heap[1] = heap[heapIndex]
				heapIndex -= 1
				x = heap[1]
				j0 = 1

				while true
					j1 = 2*j0
					j2 = j1+1
	
					if j1 > heapIndex
						break
					end
	
					if j2 > heapIndex
						k0 = j1
					else
						if heap[j1] < heap[j2]
							k0 = j1
						else
							k0 = j2
						end
					end
	
					if heap[k0] < x
						heap[j0] = heap[k0]
						j0 = k0
					else
						break
					end	
				end

				heap[j0] = x   
				u = ifts[u]  
				v = heap[1]  
				v = ifts[v] 
				
				if u == v
					edgeDiff[u] -= we
				else
					w = lcaQuery(u, v, lca_LV, lca_ET, lca_R, lca_M)
					x = fts[w]
					heapIndex += 1
					j0 = heapIndex
					while j0 > 1
						pj = Int(floor(j0/2))
						if heap[pj]> x
							heap[j0] = heap[pj]
						else
							break 
						end
						j0 = pj
					end
					heap[j0] = x
				end
			end  
	
			w = ifts[heap[1]]		
			edgeDiff[w] -= we    
			edgeTerminators[j] = w
			heapIndex = 0
		end
	end

	return (edgeDiff, edgeTerminators)
end

function computeEdgeCuts_f(H::Hypergraph, p::Pindex, lcaInstance::LCA, E::Vector{Int}, hyperedges_flag::Vector{Int})	
	w_ = H.hwts
	we = 0
	lca_R = lcaInstance.fts
	lca_ET = lcaInstance.eulerTour
	lca_LV = lcaInstance.levelVec
	lca_M = lcaInstance.rmqSparseTable
	ifts = lcaInstance.ifts
	fts = lcaInstance.fts
	n = H.n
	M = H.e
	v_ = H.hedges
	loc_ = H.eptr	
	isfixed = zeros(Int, n)
	isfixed[p.p1] .= 1
	isfixed[p.p2] .= 1
	edgeDiff = zeros(Int, n)
	edgeTerminators = zeros(Int, M)
	heap = zeros(Int, 2*max(maximum(loc_[2:end]-loc_[1:end-1]),1))  
	heapIndex = 0
	cedge = zeros(Int, 50000)

	for jj in 1:length(E) 
		j = E[jj]

        if hyperedges_flag[j] == 0
            continue
        end

		we = w_[j]
		l_edge = 0
		for kk in loc_[j]:loc_[j+1]-1
			if isfixed[v_[kk]] == 0
				l_edge += 1
				cedge[l_edge] = v_[kk]
			end    
		end
		
		if l_edge > 0  
			if l_edge == 1                     
				u = cedge[1]
				edgeTerminators[j] = u
				edgeDiff[u] += we 
			elseif l_edge == 2    
				u = cedge[1]
				v = cedge[2]
				w = lcaQuery(u, v, lca_LV, lca_ET, lca_R, lca_M)
				edgeTerminators[j] = w

				if w == u
					edgeDiff[v] += we
				elseif w == v
					edgeDiff[u] += we
				else
					edgeDiff[v] += we
					edgeDiff[u] += we
					edgeDiff[w] -= we
				end	
			else   
				for k in 1:l_edge
					u = cedge[k]
					edgeDiff[u] += we
					u = fts[u]
					x = u
					heapIndex += 1
					j0 = heapIndex

					while j0 > 1
						pj = Int(floor(j0/2))
	
						if heap[pj] > x
							heap[j0] = heap[pj]
						else
							break    
						end 
	
						j0 = pj
					end

					heap[j0] = x
				end 
	
				while heapIndex > 1
					u = heap[1]
					heap[1] = heap[heapIndex]
					heapIndex -= 1
					x = heap[1]
					j0 = 1

					while true
						j1 = 2*j0
						j2 = j1+1
	
						if j1 > heapIndex     
							break
						end
	
						if j2 > heapIndex
							k0 = j1
						else
							if heap[j1] < heap[j2]
								k0 = j1
							else
								k0 = j2
							end
						end
	
						if heap[k0] < x
							heap[j0] = heap[k0]
							j0 = k0
						else
							break
						end	
					end

					heap[j0] = x
					u = ifts[u]
					v = heap[1]
					v = ifts[v]
	
					if u == v
						edgeDiff[u] -= we 
					else
						w = lcaQuery(u, v, lca_LV, lca_ET, lca_R, lca_M)
						x = fts[w]
						heapIndex += 1
						j0 = heapIndex

						while j0 > 1
							pj = Int(floor(j0/2))
							if heap[pj]> x
								heap[j0] = heap[pj]
							else
								break    
							end
							j0 = pj
						end

						heap[j0] = x
					end
				end  

				w = ifts[heap[1]]			
				edgeTerminators[j] = w
				heapIndex = 0
			end
		end 
	end 
	return (edgeDiff, edgeTerminators)
end

function analyzeCutsOnTree(H::Hypergraph, B::Incidence, vtxWts::Vector{Int}, p::Pindex, T::SimpleWeightedGraphs.SimpleGraph, hyperedges_flag::Vector{Int})
	e = H.e
	n = H.n
	w_ = H.hwts

	forced_0 = incidentEdges(H, B, p.p1, hyperedges_flag)
	forced_1 = incidentEdges(H, B, p.p2, hyperedges_flag)
	forced_01 = intersect(forced_0, forced_1)

	setdiff!(forced_0, forced_01)
	setdiff!(forced_1, forced_01)
	forcedType = -ones(Int, e)
	forcedType[forced_0] .= 0
	forcedType[forced_1] .= 1
	forcedType[forced_01] .= 2

	(rmqSparseTable, eulerLevel, eulerVecs, ~) = lca2rmq(T, 1)

	child = eulerVecs[3]
	parent = eulerVecs[4]
	eulerTour = eulerVecs[1]
	fts = zeros(Int, n)

	for idx in 1:length(eulerTour)
        fts[eulerTour[idx]] = idx
    end

	ifts = zeros(Int, length(eulerLevel))

    for idx in 1:length(fts)
        ifts[fts[idx]] = idx
	end

	lcaInstance = LCA(rmqSparseTable, eulerLevel, child, parent, eulerTour, eulerLevel, fts, ifts)
	mobileEdges = findall(x-> x==-1, forcedType)
	(edgeDiff, edgeTerminators) = computeEdgeCuts(H, lcaInstance, mobileEdges, hyperedges_flag)
	(edgeDiff0, edgeTerminators0) = computeEdgeCuts_f(H, p, lcaInstance, forced_0, hyperedges_flag)
	(edgeDiff1, edgeTerminators1) = computeEdgeCuts_f(H, p, lcaInstance, forced_1, hyperedges_flag)

	edgeTerminators += edgeTerminators0 + edgeTerminators1
	vtxCuts = zeros(Int, n)
	edgeCuts = zeros(Int, n)
	edgeCuts0 = zeros(Int, n)
	edgeCuts1 = zeros(Int, n)
	FB0 = zeros(Int, n)
	FB1 = zeros(Int, n)

	for j in 1:length(forced_0)
		edge = forced_0[j]
		FB0[edgeTerminators[edge]] += w_[edge]
	end
	
	for j in 1:length(forced_1)
		edge = forced_1[j]
		FB1[edgeTerminators[edge]] += w_[edge]
	end    

	for j in 1:length(child)-1 
		nodec = child[j]
		pnodec = parent[nodec]
		edgeCuts[nodec] += edgeDiff[nodec]
		edgeCuts[pnodec] += edgeCuts[nodec]	
		edgeCuts0[nodec] +=  edgeDiff0[nodec]
		edgeCuts0[pnodec] += edgeCuts0[nodec]
		edgeCuts1[nodec] += edgeDiff1[nodec]
		edgeCuts1[pnodec] += edgeCuts1[nodec]
		vtxCuts[nodec] += vtxWts[nodec]
		vtxCuts[pnodec] += vtxCuts[nodec]
		FB0[pnodec] += FB0[nodec]
		FB1[pnodec] += FB1[nodec]
	end

	return CutProfile(vtxCuts, edgeCuts, edgeDiff, parent, edgeTerminators, p, 
				   forcedType, forced_0, forced_1, forced_01, 
				   FB0, FB1, edgeCuts0, edgeCuts1)
end
