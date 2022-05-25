function IsolateIslands(H::Hypergraph)
    cid = 0
    nidx = 0
    ptr = 0
    ptr1 = 0
    ptr2 = 0
    n = H.n
    e = H.e
    hedges = H.hedges
    eptr = H.eptr
    vwts = H.vwts
    hwts = H.hwts
    vwts_new = similar(vwts)
    hwts_new = similar(hwts)
    k =  length(hedges)
    originalIndices = zeros(Int, n)
    unusedIndices = zeros(Int, n)
    newIndices = zeros(Int, n)
    hedges_new = zeros(Int, k)
    eptr_new = zeros(Int, length(eptr))
    (cID, csizes) = hypergraphCC(H, Vector{Int}(), false)
    id = findmax(csizes)[2] 

    for i in 1:length(cID)
        if cID[i] == id
            ptr1 += 1
            originalIndices[ptr1] = i
            newIndices[i] = ptr1
            vwts_new[ptr1] = vwts[i]
        else
            ptr2 += 1
            unusedIndices[ptr2] = i
        end 
    end

    originalIndices = originalIndices[1:ptr1]
    vwts_new = vwts_new[1:ptr1]
    unusedIndices = unusedIndices[1:ptr2]

    for i in 1:k
        hedge = hedges[i]
        cid = cID[hedge]
        if cid   == id
            ptr += 1
            nidx = newIndices[hedge]
            hedges_new[ptr] = nidx
        end
    end

    hedges_new = Vector(hedges_new[1:ptr])
    ptr = 0
    x = Int(1)
    ne = Int(0)
    size = Int(0)
    nv = maximum(newIndices)    

    for k in 1:e
        loc = eptr[k]
        if cID[hedges[loc]] == id
            ne += 1
            ptr += 1
            size = eptr[k+1]-eptr[k]
            eptr_new[ptr] = x
            hwts_new[ptr] = hwts[k]
            x += size
        end
    end
    
    ptr += 1    
    eptr_new[ptr] = x
    eptr_new = Vector(eptr_new[1:ptr])  
    hwts = hwts_new[1:ptr-1]
    
    return (Hypergraph(nv, ne, hedges_new, eptr_new, vwts_new, hwts_new), originalIndices, newIndices, unusedIndices)
end