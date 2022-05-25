#=function ReadHmetisBenchmark(fname::String)
    i = 0
    j = 0
    thr = 3
    hE = 0
    n = 0
    l = 1
    eptr = Int[]
    hedges = Int[]
    vwts = Int[]
    f = open(fname)
    vweight = false

    for ln in eachline(f)
        i += 1

        if i > 1 && i <= thr
            matchline = match(r"([0-9\s]+)", ln)

            if matchline != nothing
                hdge = [parse(Int, vtxstr) for vtxstr in split(matchline.captures[1])]
                append!(hedges, hdge)
                l += length(hdge)
                eptr[i] = l
            end

        elseif i > thr
            vweight = true
            j += 1
            vwts[j] = parse(Int, ln)
        else
            match_first_line = match(r"(\d+)\s+(\d+)", ln)
            hE = parse(Int, match_first_line.captures[1])
            n = parse(Int, match_first_line.captures[2])
            resize!(eptr, hE + 1)
            resize!(vwts, n)
            eptr[1] = 1
            thr = hE+1
        end
    end

    close(f)

    if vweight == false
        vwts = ones(Int, n)
    end

    return (hedges, eptr, vwts, n, hE, length(hedges))
end=#

function ReadHypergraphFile(fname::String)
    i = 0
    j = 0
    num_hyperedges = 0
    num_vertices = 0
    wt_type = 0
    l = 1
    eptr = Int[]
    hedges = Int[]
    hyperedge_wts = Int[]
    vertex_wts = Int[]
    f = open(fname, "r")

    for ln in eachline(f)
        if i == 0
            line_vector = split(ln)
            num_vertices = parse(Int, line_vector[2])
            num_hyperedges = parse(Int, line_vector[1])
            eptr = zeros(Int, num_hyperedges+1)
            eptr[1] = 1
            hyperedge_wts = zeros(Int, num_hyperedges)
            vertex_wts = ones(Int, num_vertices)

            if length(line_vector) > 2
                wt_type = parse(Int, line_vector[3])
            end
        elseif i <= num_hyperedges
            line_vector = split(ln)
            edge_vector = parse.(Int, line_vector)

            if wt_type == 10 || wt_type == 0
                hyperedge_wts[i] = 1
                append!(hedges, edge_vector)
                l += length(edge_vector)
                eptr[i+1] = l
            else
                hyperedge_wts[i] = edge_vector[1]
                append!(hedges, edge_vector[2:end])
                l += length(edge_vector)-1
                eptr[i+1] = l
            end
        else
            vertex_wts[i-num_hyperedges] = parse(Int, ln)
        end

        i += 1
    end

    close(f)

    return (hedges, eptr, vertex_wts, hyperedge_wts, num_vertices, num_hyperedges, length(hedges))
end

function ReadHypergraphFixedFile(fname::String, fixed_vertices::Vector{Int})
    i = 0
    f = open(fname, "r")

    for ln in eachline(f)
        i += 1
        fixed_vertices[i] = parse(Int, ln)
    end
end