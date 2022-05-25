function GenerateHyperedgesHash(H::Hypergraph)
    hyperedges_hash = Dict{Int, Int}()

    for i in 1:H.e
        start_idx = H.eptr[i]
        end_idx = H.eptr[i+1]-1
        hash_score = 0

        for j in start_idx:end_idx
            v = H.hedges[j]
            hash_score += v*v
        end

        push!(hyperedges_hash, i => hash_score)
    end

    return hyperedges_hash
end
