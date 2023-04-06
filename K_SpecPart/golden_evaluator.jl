function golden_evaluator(hypergraph::__hypergraph__,
                        num_parts::Int,
                        partition::Vector{Int})
    cutsize = 0
    for i in 1:hypergraph.num_hyperedges
        first_valid_entry = hypergraph.eptr[i]
        first_invalid_entry = hypergraph.eptr[i+1]
        for j in first_valid_entry+1 : first_invalid_entry-1
            v = hypergraph.eind[j]
            if (partition[hypergraph.eind[j]] != partition[hypergraph.eind[j-1]])
                cutsize += hypergraph.hwts[i]
                break
            end
        end
    end
    balance = zeros(Int, num_parts)
    for i in 1:hypergraph.num_vertices
        p = partition[i]
        balance[p+1] += hypergraph.vwts[i]
    end
    return (cutsize, balance)
end