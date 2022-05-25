function GenerateSpectralCommunities(bins::Matrix{Int}, hypergraph_c::Hypergraph_C, incidence_struct::Incidence)
    hypergraph = hypergraph_c.HG
    vwts = hypergraph.vwts
    fixed = hypergraph_c.fixed_part
    cc, cs = ClusterHG(bins, hypergraph, incidence_struct)
    vwts_cc = ContractVtxWts(vwts, cc)
    n_cc, e_cc, hedges_cc, eptr_cc, hwts_cc = ContractHyperGraph(hypergraph, cc)
    hypergraph_cc = Hypergraph(n_cc, e_cc, hedges_cc, eptr_cc, vwts_cc, hwts_cc)
    incidence_struct_cc = HypergraphToIncidence(hypergraph_cc)
    fixed_cc = -ones(Int, hypergraph_cc.n)

    for i in 1:length(fixed)
        if fixed[i] > -1
            fixed_cc[cc[i]] = fixed[i]
        end
    end

    return (hypergraph_cc, incidence_struct_cc, fixed_cc, vwts_cc, cc, cs)
end

function hypergraphCC(H::Hypergraph, excludedEdges::Vector{Int}, mode::Bool)
    eptr = H.eptr
    hedges = H.hedges
    e = H.e
    n = H.n
    id = Vector{Int}(1:n)
    sz = ones(Int, n)

    includedEdges = setdiff(Vector{Int}(1:e), excludedEdges)

    for k in 1:length(includedEdges)
        kk = includedEdges[k]
        l_edge = eptr[kk+1]-eptr[kk]

        for j in 0:l_edge-2
            u = hedges[eptr[kk]+j]
            v = hedges[eptr[kk]+j+1]

            while u != id[u]
                id[u] = id[id[u]]
                u = id[u]
            end

            while v != id[v]
                id[v] = id[id[v]]
                v = id[v]
            end

            if u != v
                if sz[u] < sz[v]
                    id[u] = v
                    sz[v] += sz[u]
                else
                    id[v] = u
                    sz[u] += sz[v]
                end
            end
        end
    end

    cID = zeros(Int, n)
    cs = zeros(Int, n)

    for k in 1:n
        u = k
        
        while u != id[u]
            u = id[u]
        end

        cID[k] = u
        cs[u] += 1
    end

    ndx = findall(x-> x > 0, cs)
    ord = zeros(Int, maximum(ndx))
    ord[ndx] = Vector{Int}(1:length(ndx))

    for k in 1:n
        cID[k] = ord[cID[k]]
    end

    cs = cs[ndx]

    return (cID, cs)
end

function standardizeClustering(cID::Vector{Int}, cs::Vector{Int})
	cID .+= 2
    cID[1] = 1
    cID[2] = 2
    B = zeros(Int, length(cID))

    @inbounds for j in 1:length(cID)
        B[cID[j]] = 1
    end

    X = findall(!iszero, B)
    M = zeros(Int, maximum(X))
    M[X] = 1:length(X)
    cID = M[cID]
    cs[3:end] = cs[1:end-2]
    cs[1:2] .= 1

    return (cID, cs)
end

#=
function ContractHyperGraph(H::Hypergraph, cID::Vector{Int})
    c_ptr = 1
    n = maximum(cID)

    for i in 1:H.eptr
        active = false
        start_idx = H.eptr[i]
        end_idx = H.eptr[i+1]
        vtxs = H.hedges[start_idx:end_idx-1]
        base_c = cId[H.vtxs[1]]

        for j in start_idx+1:end_idx
            cid = vtxs[j]
            if cid != base_c
                break
            else
                active = true
            end
        end
end
=#

function ContractHyperGraph(H::Hypergraph, cID::Vector{Int})
    hedges = H.hedges
    eptr = H.eptr
    e = H.e

    e_ = zeros(Int, length(hedges)+1)
    e_[eptr] .= 1
    e_ = cumsum(e_)
    e_ = e_[1:end-1]

    v_ = cID[hedges]
    E = sparse(v_, e_, ones(Int, length(e_)))

    (v_, e_, ~) = findnz(E)
    E = sparse(v_, e_, ones(Int, length(e_)))

    E_sum = sum(E, dims=1)[1, :]

    ndx = findall(x-> x>1, E_sum)
    E = E[:,ndx] 

    w_ = H.hwts
    w_ = w_[ndx]

    (n, ~) = size(E)
    r = rand(1, n)
    h = (r*E)[1, :]
    prm = sortperm(h)
    hs = h[prm] 

    ind = findall(hs[1:end-1] .!= hs[2:end])
    push!(ind, length(hs))

    w_ = w_[prm]
    w1 = cumsum(w_)
    w1[ind[2:end]] -= w1[ind[1:end-1]]

    ind = prm[ind]
    w1[prm] = w1
    E = E[:,ind]      
    (v_, e_, ~) = findnz(E)  
    w_ = w1[ind]

    loc_ = findall(e_[1:end-1] .!= e_[2:end])
    push!(loc_, length(e_))
    push!(loc_, 0)
    loc_[2:end] = loc_[1:end-1] .+ 1
    loc_[1] = 1

    return (maximum(v_), length(loc_)-1, v_, loc_, w_)
end

function ContractVtxWts(vwts::Vector{Int}, cc::Vector{Int})
    cmax = maximum(cc)
    vwts_c = zeros(Int, cmax)

    for i in 1:length(cc)
        vwts_c[cc[i]] += vwts[i]
    end

    return vwts_c
end

function ClusterHG(bins::Matrix{Int}, H::Hypergraph, B::Incidence)
    (m, ~) = size(bins)
    union_cut = Int[]
    intersect_cut = Vector{Int}(1:H.e)
    hyperedges = Vector{Int}(1:H.e)

    for i in 1:m
        bins_i = bins[i, :]
        (cut_edges_mrk, ~, ~, ~, ~) = cutProfile(H, B, bins_i)
        cut_edges_i = findall(!iszero, cut_edges_mrk)
        union!(union_cut, cut_edges_i)
        intersect!(intersect_cut, cut_edges_i)
    end

    (cc, cs) = hypergraphCC(H, union_cut, false)

    return cc, cs
end
