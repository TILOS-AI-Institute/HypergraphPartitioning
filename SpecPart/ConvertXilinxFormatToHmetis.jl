function ConvertXilinxToHmetis(hypergraph_file::String)
    line_idx = 0
    l = 1
    netlist_file = hypergraph_file * "/HMGraph_0"
    nodelist_file = hypergraph_file * "/HMNodeInfo_0"
    num_vertices = 0
    num_hyperedges = 0
    hyperedge_wts = Vector{Int}()
    vertex_wts = Vector{Int}()
    hedges = Vector{Int}()
    eptr = Vector{Int}()
    fixed_vtxs = Vector{Int}()

    f = open(netlist_file, "r")

    for ln in eachline(f)
        if line_idx == 0
            m_first_line = match(r"NumEdges : (\d+) NumNodes : (\d+)", ln)
            num_hyperedges = parse(Int, m_first_line.captures[1])
            num_vertices = parse(Int, m_first_line.captures[2])
            hyperedge_wts = zeros(Int, num_hyperedges)
            vertex_wts = zeros(Int, num_vertices)
            fixed_vtxs = zeros(Int, num_vertices)
            eptr = zeros(Int, num_hyperedges+1)
            eptr[1] = 1
        elseif line_idx == 1
            line_idx += 1
            continue
        else
            m_edge_list = match(r"(\d+)\s+(\d+)\s+:\s+([\d\s]+)", ln)
            hyperedge_id = parse(Int, m_edge_list.captures[1])
            hyperedge_wt = parse(Int, m_edge_list.captures[2])
            hyperedge_wts[hyperedge_id+1] = hyperedge_wt
            vertex_vector = [parse(Int, vtxstr) for vtxstr in split(m_edge_list.captures[3])] .+ 1
            append!(hedges, vertex_vector)
            l += length(vertex_vector)
            eptr[hyperedge_id+2] = l
        end

        line_idx += 1
    end

    close(f)

    line_idx = 0

    f = open(nodelist_file, "r")

    for ln in eachline(f)
        if line_idx <= 1
            line_idx += 1
            continue
        end

        m_vertex_info = match(r"(\d+)\s+(\d+)\s+(-?\d+)\s+(\d+)\s+:\s+(\d+)", ln)
        vertex_id = parse(Int, m_vertex_info.captures[1])
        fixed_info = parse(Int, m_vertex_info.captures[3])
        fixed_vtxs[vertex_id+1] = fixed_info
        vtx_wt = parse(Int, m_vertex_info.captures[5])
        vertex_wts[vertex_id+1] = vtx_wt
        line_idx += 1
    end

    close(f)

    return (Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_wts, hyperedge_wts), fixed_vtxs)
end