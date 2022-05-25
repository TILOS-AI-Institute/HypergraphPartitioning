function PlotEmbedding(hypergraph_file::String, h_c::Hypergraph_C; eigenvecs::Int = 2, expander_cycles::Int = 3, eigen_iters::Int = 20, bestsolns_on_tree::Int = 10)
    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    (adj_mat, ~) = HypergraphToGraph(hypergraph, expander_cycles)
    X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)
    hmetis(hypergraph_file, 2, 5, 10, 1, 1, 0, 1, 24)
    golden_vertices = zeros(Int, hypergraph.n)
    hypergraph_part = hypergraph_file*".part.2"
    f = open(hypergraph_part, "r")
    i = 0

    for ln in eachline(f)
        i += 1
        golden_vertices[i] = parse(Int, ln)
    end

    f = open(hypergraph_file*".embedding.dat", "w")
    df = DataFrame(:Eigen1 => X[:, 1], :Eigen2 => X[:, 2], :Partition => golden_vertices)
    CSV.write(hypergraph_file*".embedding.dat", df)
    close(f)
end