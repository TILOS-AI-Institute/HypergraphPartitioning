function ExportHypergraph(hypergraph::Hypergraph, fname::String, fixed_vtxs::Vector{Int}; idx_flag::Int = 1, wt_flag::Int = 0)
    n = hypergraph.n
    e = hypergraph.e
    hedges = hypergraph.hedges
    eptr = hypergraph.eptr
    hwts = hypergraph.hwts
    vwts = hypergraph.vwts

    f = open(fname, "w")

    println(f, e, " ", n, " 11")

    for i in 1:e
        start_idx = eptr[i]
        end_idx = eptr[i+1]-1
        print(f, hwts[i])

        if idx_flag == 1
            for j in start_idx:end_idx
                print(f, " ", hedges[j])
            end
            print(f, "\n")
        else
            for j in start_idx:end_idx
                print(f, " ", hedges[j]-1)
            end
            print(f, "\n")
        end
    end

    for i in 1:n
        if wt_flag == 0
            println(f, vwts[i])
        else
            println(f, vwts[i]+1)
        end
    end

    close(f)

    if maximum(fixed_vtxs) > -1
        f = open(fname*".fixed", "w")

        for i in 1:n
            println(f, fixed_vtxs[i])
        end

        close(f)
    end
end