import Graphs

function BuildMetisGraph(tree::SimpleWeightedGraph, opts::String)
    time_stamp = string(now())
    fname = "metis_graph" * opts * "." * time_stamp * ".gr"
    f = open(fname, "w")
    wts = tree.weights
    println(f, SimpleWeightedGraphs.nv(tree), " ", SimpleWeightedGraphs.ne(tree), " 001" )

    for i in 1:SimpleWeightedGraphs.nv(tree)
        #nbrs = SimpleWeightedGraphs.neighbors(tree, i)
        nbrs = Graphs.neighbors(tree, i)

        for j in 1:length(nbrs)
            nbr_vtx = nbrs[j]
            wt = Int(wts[i, nbr_vtx])
            #wt = Int(round(tree.weights[i, nbr_vtx] * 1e07))
            #wt = wt == 0 ? 1 : wt
            print(f, nbr_vtx, " ", wt, " ")
        end

        print(f, "\n")
    end

    close(f)

    return fname
end
