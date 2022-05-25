function MSDR(emst::SimpleWeightedGraph)
    i = 0
    eps = 0.0001

    edge_list = zeros(Int, nv(emst), 2)
    edge_weight_list = zeros(ne(emst))

    for e in edges(emst)
        i += 1
        edge_list[i, :] = [e.src, e.dst]
        edge_weight_list[i] = e.weight
    end

    edge_stdev = std(edge_weight_list)
    sub_trees = deepcopy(emst)
    i = 0

    while done == false
        i += 1
        temp = edge_stdev

        for e in edges(sub_trees)



end