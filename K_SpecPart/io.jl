using DelimitedFiles

function read_hypergraph_file(file_name::String, 
                            fixed_file_name::String = "")
    num_vertices = 0
    num_hyperedges = 0
    line_count = 0
    eptr_count = 1
    wt_type = 0
    eptr = Int[]
    eind = Int[]
    vwts = Int[]
    hwts = Int[]
    fixed_vertices = Int[]
    f = open(file_name, "r")

    for ln in eachline(f)
        if line_count == 0
            line_vector = split(ln)
            num_vertices = parse(Int, line_vector[2])
            num_hyperedges = parse(Int, line_vector[1])
            eptr = zeros(Int, num_hyperedges + 1)
            eptr[1] = 1
            hwts = ones(Int, num_hyperedges)
            vwts = ones(Int, num_vertices)
            fixed_vertices = -ones(Int, num_vertices)
            if length(line_vector) > 2
                wt_type = parse(Int, line_vector[3])
            end
        elseif line_count <= num_hyperedges
            line_vector = split(ln)
            hyperedge = parse.(Int, line_vector)
            if wt_type == 10 || wt_type == 0
                append!(eind, hyperedge)
                eptr_count += length(hyperedge)
                eptr[line_count+1] = eptr_count
            else
                hwts[line_count] = hyperedge[1]
                append!(eind, hyperedge[2:end])
                eptr_count += length(hyperedge) - 1
                eptr[line_count+1] = eptr_count
            end
        else
            vwts[line_count-num_hyperedges] = parse(Int, ln)
        end

        line_count += 1
    end

    close(f)

    if fixed_file_name != ""
        read_fixed_vertices_file(fixed_file_name, fixed_vertices)
    end
    return build_hypergraph(num_vertices, 
                        num_hyperedges, 
                        eptr, 
                        eind, 
                        fixed_vertices, 
                        vwts, 
                        hwts)
end

function read_fixed_vertices_file(file_name::String, 
                                fixed_vertices::Vector{Int}) 
    line_count = 0
    f = open(file_name, "r")
    
    for ln in eachline(f)
        line_count += 1
        fixed_vertices[line_count] = parse(Int, ln)
    end
end

function read_hint_file(file_name::String)
    line_count = 0
    p = readdlm(file_name, Int)
    partition = p[:,1]
    return partition
    #=f = open(file_name, "r")
    for ln in eachline(f)
        line_count += 1
        partition[line_count] = parse(Int, ln)
    end=#
end