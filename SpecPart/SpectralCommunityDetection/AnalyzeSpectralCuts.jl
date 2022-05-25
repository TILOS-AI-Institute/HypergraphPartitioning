function CompareCutHyperedges(spectral_cut::Vector{Int}, hmetis_cut::Vector{Int}, hypergraph::Hypergraph, incidence_struct::Incidence)
    (spectral_cut_edges_mrk, ~, ~, ~, ~) = cutProfile(hypergraph, incidence_struct, spectral_cut)
    spectral_cut_hyperedges = findall(!iszero, spectral_cut_edges_mrk)

    (hmetis_cut_edges_mrk, ~, ~, ~, ~) = cutProfile(hypergraph, incidence_struct, hmetis_cut)
    hmetis_cut_hyperedges = findall(!iszero, hmetis_cut_edges_mrk)

    hyperedges_overlap = intersect(spectral_cut_hyperedges, hmetis_cut_hyperedges)

    @info "OVERLAPPING HYPEREDGES: $(length(hyperedges_overlap))"
end

function ReadHmetisCut(fname::String, hypergraph::Hypergraph)
    hmetis(fname, 2, 5, 10, 1, 1, 0, 1, 24)

    f = open(fname*".part.2", "r")
    i = 0
    pvec = zeros(Int, hypergraph.n)

    for ln in eachline(f)
        i += 1
        pvec[i] = parse(Int, ln)
    end
    
    close(f)

    return FindCutSize(pvec, hypergraph), pvec
end