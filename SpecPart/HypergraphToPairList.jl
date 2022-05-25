function GenerateHypergraphPairList(hypergraph::Hypergraph)
    hyperedges_pair_list = zeros(Int, hypergraph.e, 2)

    for i in 1:hypergraph.e
        start_idx = hypergraph.eptr[i]
        end_idx = hypergraph.eptr[i+1]
        hyperedges_pair_list[i, 1] = start_idx
        hyperedges_pair_list[i, 2] = end_idx-1
    end

    return hyperedges_pair_list
end
