function GenerateHyperstretch(H::Hypergraph, X::Array{Float64})
    hyperstretch_hyperedges = zeros(H.e)
    (~, m) = size(X)

    for i in 1:H.e
        start_idx = H.eptr[i]
        end_idx = H.eptr[i+1]-1
        vtxs = H.hedges[start_idx:end_idx]
        vtxs_embedding = X[vtxs, :]
        average_embedding = mean(vtxs_embedding, dims=1)
        hyperstretch_vector = abs.(average_embedding .- vtxs_embedding)
        hyperstretch = sum(hyperstretch_vector, dims=1)
        total_stretch = 0.0

        for j in 1:m
            total_stretch += hyperstretch[j] * hyperstretch[j]
        end

        hyperstretch_hyperedges[i] = sqrt(total_stretch)
    end

    return hyperstretch_hyperedges
end