function cutProfile(H::Hypergraph, B::Incidence, bins::Vector{Int})
    n = H.n
    e = H.e
    eptr = H.eptr
    w_ = H.hwts
    hedges = H.hedges
    nIncidentEdgesCut = zeros(Int, n)
    cutEdges = zeros(Int, e)
    cutStatistics = zeros(Int, e)
    edgePart = zeros(Int, e)
    fp = zeros(Int, n)
    S1 = zeros(Int, n)
    nodesInCutOnly = ones(Int, n)
    b_hedges = B.hedges
    b_eptr = B.eptr

    for j in 1:e
        cutFlag = false
        
        for k in eptr[j]:eptr[j+1]-2
            if bins[hedges[k]] != bins[hedges[k+1]] 
                cutFlag = true
                break
            end
        end

        if cutFlag == true    
            edge = hedges[eptr[j]:(eptr[j+1]-1)]
            hye_parts = bins[edge]
            cutEdges[j] = 1
            edgePart[j] = -1
            s = sum(hye_parts)
        
            if s == length(edge)-1
                id = findall(iszero, hye_parts)
                fp[edge[id]] .+= w_[j]
            end
            
            if s == 1
                id = findall(x->x==1, hye_parts)
                fp[edge[id]] .+= w_[j]
            end
        else
            edgePart[j] = bins[hedges[eptr[j]]] 
        end
    end
    
    for j in 1:n
        for k in b_eptr[j]:b_eptr[j+1]-1
            if cutEdges[b_hedges[k]] == 0
                nodesInCutOnly[j] = 0
                break
            end
        end
    end

    ce = findall(!iszero, cutEdges)
    ww = 0

    for j in 1:length(ce)
        edge = ce[j]
        ww = w_[edge]
        
        for k in eptr[edge]:(eptr[edge+1]-1)
            nIncidentEdgesCut[hedges[k]] = nIncidentEdgesCut[hedges[k]] + ww
        end
    end

    return (cutEdges, nodesInCutOnly, nIncidentEdgesCut, fp, edgePart)
end

function findCutsize(bins::Vector{Int}, H::Hypergraph, B::Incidence)
    cutsize = 0
    (cutEdges, ~, ~, ~, ~) = cutProfile(H, B, bins)
    cutsize = sum(H.hwts[findall(!iszero, cutEdges)])

    return cutsize
end
