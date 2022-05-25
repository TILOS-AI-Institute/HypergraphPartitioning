function StructureDetection(h_c::Hypergraph_C, incidence_struct::Incidence, ubfac::Int, eigs::Int; community_opts::Int=1)
    community = zeros(Int, h_c.HG.n)
   
    if community_opts == 1
        community = SpectralCommunities(h_c, incidence_struct, ubfac, eigenvecs = eigs)
    end

    return community
end
