function findLabels(clusters::AbstractArray, n::Int)
	label = zeros(Int, n)

	@inbounds for i in 1:length(clusters)
		for j in 1:length(clusters[i])
			label[clusters[i][j]] = i
		end
	end
    
	return label
end

function incidentEdges(H::Hypergraph, B::Incidence, S::Vector{Int})
    hedges = B.hedges
	eptr = B.eptr
	e = H.e
	T = 0
	
    @inbounds for j in 1:length(S)
        k = S[j]
        T += eptr[k+1]-eptr[k]
    end    

    E = zeros(Int, T)
	r = 1
	
    @inbounds for j in 1:length(S)
        k = S[j]
        L = eptr[k+1]-eptr[k]
        E[r:r+L-1] = hedges[eptr[k]:(eptr[k+1]-1)]
        r += L;
    end      
    
	unique!(E)
	return E
end

function incidentNodes(H::Hypergraph, E::Vector{Int})
    n = H.n
    eptr = H.eptr
    hedges = H.hedges
    
	T = 0
	
    @inbounds for j in 1:length(E)
        k = E[j]
        T += eptr[k+1]-eptr[k]
    end

    S = zeros(Int, T)
	r = 1
	
    @inbounds for j in 1:length(E)
        k = E[j]
        L = eptr[k+1]-eptr[k]
        S[r:r+L-1] = hedges[eptr[k]:(eptr[k+1]-1)]
        r += L 
    end      
    
    unique!(S)
	return S
end

function computeEdgeCuts(H::Hypergraph, lcaInstance::LCA, E::Vector{Int})
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

function computeEdgeCuts_f(H::Hypergraph, p::Pindex, lcaInstance::LCA, E::Vector{Int})	
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

	@inbounds for jj in 1:length(E) 
		j = E[jj]
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

function analyzeCutsOnTree(H::Hypergraph, B::Incidence, vtxWts::Vector{Int}, p::Pindex, T::AbstractGraph)
	e = H.e
	n = H.n
	w_ = H.hwts

	forced_0 = incidentEdges(H, B, p.p1)
	forced_1 = incidentEdges(H, B, p.p2)
	forced_01 = intersect(forced_0, forced_1)

	setdiff!(forced_0, forced_01)
	setdiff!(forced_1, forced_01)
	forcedType = -ones(Int, e)
	@inbounds forcedType[forced_0] .= 0
	@inbounds forcedType[forced_1] .= 1
	@inbounds forcedType[forced_01] .= 2

	(rmqSparseTable, eulerLevel, eulerVecs, ~) = lca2rmq(T, 1)

	child = eulerVecs[3]
	parent = eulerVecs[4]
	eulerTour = eulerVecs[1]
	fts = zeros(Int, n)

	@inbounds @simd for idx in 1:length(eulerTour)
        fts[eulerTour[idx]] = idx
    end

	ifts = zeros(Int, length(eulerLevel))

    @inbounds @simd for idx in 1:length(fts)
        ifts[fts[idx]] = idx
	end

	lcaInstance = LCA(rmqSparseTable, eulerLevel, child, parent, eulerTour, eulerLevel, fts, ifts)
	mobileEdges = findall(x-> x==-1, forcedType)
	(edgeDiff, edgeTerminators) = computeEdgeCuts(H, lcaInstance, mobileEdges)
	(edgeDiff0, edgeTerminators0) = computeEdgeCuts_f(H, p, lcaInstance, forced_0)
	(edgeDiff1, edgeTerminators1) = computeEdgeCuts_f(H, p, lcaInstance, forced_1)

	edgeTerminators += edgeTerminators0 + edgeTerminators1
	vtxCuts = zeros(Int, n)
	edgeCuts = zeros(Int, n)
	edgeCuts0 = zeros(Int, n)
	edgeCuts1 = zeros(Int, n)
	FB0 = zeros(Int, n)
	FB1 = zeros(Int, n)

	@inbounds for j in 1:length(forced_0)
		edge = forced_0[j]
		FB0[edgeTerminators[edge]] += w_[edge]
	end
	
	@inbounds for j in 1:length(forced_1)
		edge = forced_1[j]
		FB1[edgeTerminators[edge]] += w_[edge]
	end    

	@inbounds for j in 1:length(child)-1 
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

function BestCutsOnTree(T::AbstractGraph, H::Hypergraph, B::Incidence, p::Pindex, capacity::Vector{Int}, k::Int)
    vtxWts = H.vwts
    binMatrix = zeros(Int, k, H.n)
	cutInfo = analyzeCutsOnTree(H, B, vtxWts, p, T)
	w_ = H.hwts
	n = H.n
	e = H.e
    vtxCuts = cutInfo.vtxCuts
	edgeCuts = cutInfo.edgeCuts
	pred = cutInfo.pred
	pindex = cutInfo.pindex
	p1 = pindex.p1
	p2 = pindex.p2
	forced_0 = cutInfo.forced_0
	forced_1 = cutInfo.forced_1
	forced_01 = cutInfo.forced_01
	FB0 = cutInfo.FB0
	FB1 = cutInfo.FB1
	edgeCuts0 = cutInfo.edgeCuts0
	edgeCuts1 = cutInfo.edgeCuts1
	nforced_0 = sum(w_[forced_0]) 
	nforced_1 = sum(w_[forced_1]) 
	nforced_01 = sum(w_[forced_01]) 
    s = sum(vtxWts)
	edgeCuts[1] = e
	vtxCuts[1] = 0
	exc0 = zeros(Int, n)
	exc1 = zeros(Int, n)
    areaPart = zeros(Int, 2)
	cutCost0 = zeros(Float64, n)
	cutCost1 = zeros(Float64, n)
	ratioCost = zeros(Float64, n)
	cutCost = zeros(Float64, n)
	areaCost = zeros(Float64, n)
	tCost_v = zeros(Float64, n)
	polarity = zeros(Int, n)
	statusFlag = zeros(Int, n)
	areaUtil0 = zeros(Int, n)
	areaUtil1 = zeros(Int, n)

	for i in 1:n
		exc0[i] = edgeCuts[i] + nforced_0 - FB0[i] + edgeCuts1[i] + nforced_01
		exc1[i] = edgeCuts[i] + nforced_1 - FB1[i] + edgeCuts0[i] + nforced_01
	end

	exc0[1] = exc1[1] = e

	for i in 1:length(exc0)
		cutCost0[i] = exc0[i]
		cutCost1[i] = exc1[i]
		(cutCost[i], pol) = findmin([cutCost0[i], cutCost1[i]])
		polarity[i] = pol-1

		if pol == 0
			areaUtil0[i] = vtxCuts[i]	
			areaUtil1[i] = s - vtxCuts[i]
            
		else
			areaUtil1[i] = vtxCuts[i]
			areaUtil0[i] = s - vtxCuts[i]
		end

        if areaUtil0[i] > capacity[1] || areaUtil1[i] > capacity[2]
            areaCost[i] = 1e09
        end

		ratioCost[i] = cutCost[i]/(areaUtil0[i] + areaUtil1[i])
		tCost_v[i] = cutCost[i] + areaCost[i]
	end
	
	cut_idx = sortperm(tCost_v)

	if tCost_v[cut_idx[1]] >= 1e09
		cut_idx = sortperm(ratioCost)
		cut_point = cut_idx[1]
		rem_edge!(T, cut_point, pred[cut_point])
		overflow = areaUtil0[cut_point] > capacity[1] || areaUtil1[cut_point] > capacity[2]
		i = 2

		while overflow == true
			cut_point = cut_idx[i]
			rem_edge!(T, cut_point, pred[cut_point])
			comps = connected_components(T)

			for j in 1:length(comps)
				vtxs = comps[j]
				total_wt = sum(H.vwts[vtxs])
				if total_wt > capacity[1]
					overflow = true
					break
				else
					overflow = false
				end
			end

			i += 1
		end

		cc = findLabels(connected_components(T), H.n)

		vwts_cc = ContractVtxWts(H.vwts, cc)
    	n_cc, e_cc, hedges_cc, eptr_cc, hwts_cc = ContractHyperGraph(H, cc)
    	H_cc = Hypergraph(n_cc, e_cc, hedges_cc, eptr_cc, vwts_cc, hwts_cc)

		println(H_cc)

		fixed_cc = -ones(Int, H.n)
		partition_vector = PartitionILP(H_cc, fixed_cc, capacity)
		
		cutsize, part_area = FindCutSize(partition_vector, H_cc)

		@info "[TREE CUT] BEST CUT RECORDED ON TREE: $cutsize WITH AREA SPLIT: $(part_area[1]) and $(part_area[2])"

		return partition_vector'
	end

	@info "[TREE CUT] BEST CUT RECORDED ON TREE: $(Int(cutCost[cut_idx[1]])) WITH AREA SPLIT: $(areaUtil0[cut_idx[1]]) and $(areaUtil1[cut_idx[1]])"

	
	#=for i in 1:n
		rem_edge!(T, i, pred[i])
		bins_o = LightGraphs.connected_components(T)
		bins = findLabels(bins_o, n) .- 1
		cutsize, part_area = FindCutSize(bins, H)
		println("[ VTX ", i, " PRED ", pred[i], "] :: Cut size calculated: ", cutsize, " with split: ", part_area, " :: Tree cut calculated: ", cutCost[i], " with split: ", areaUtil0[i], ",", areaUtil1[i])
		add_edge!(T, i, pred[i])
	end=#


	for i in 1:k
		cutpoint = cut_idx[i]
		tCost = tCost_v[cutpoint]

		if tCost >= 1e06
			binMatrix = binMatrix[1:i,:]
			break
		end

		rem_edge!(T, cutpoint, pred[cutpoint])
		bins_o = LightGraphs.connected_components(T)
		bins = findLabels(bins_o, n)

        
		@inbounds for i in 1:length(p1)
			statusFlag[p1[i]] = 1
		end

		@inbounds for i in 1:length(p2)
			statusFlag[p2[i]] = 2
		end

		cutpointflag = bins[cutpoint] == 2 ? 1 : 0

		@inbounds for i in 1:length(bins)
			if statusFlag[i] == 1
				bins[i] = 0
			elseif statusFlag[i] == 2
				bins[i] = 1
			else
				if polarity[cutpoint] == 0
					if cutpointflag == 1
						bins[i] = bins[i]%2
					else
						bins[i] -= 1
					end
				else
					if cutpointflag == 0
						bins[i] = bins[i]%2
					else
						bins[i] -= 1
					end
				end
			end
		end
        

		binMatrix[1, :] = bins
		add_edge!(T, cutpoint, pred[cutpoint])
	end

	return binMatrix

end
