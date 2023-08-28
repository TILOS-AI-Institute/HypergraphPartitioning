function build_hypergraph(num_vertices::Int,
                        num_hyperedges::Int,
                        eptr::Vector{Int},
                        eind::Vector{Int},
                        fixed_vertices::Vector{Int},
                        vwts::Vector{Int},
                        hwts::Vector{Int})
    vertices_list = [Vector{Int}() for i in 1:num_vertices]
    for i in 1:num_hyperedges
        first_valid_entry = eptr[i]
        first_invalid_entry = eptr[i+1]
        for j in first_valid_entry:first_invalid_entry-1
            v = eind[j]
            push!(vertices_list[v], i)
        end
    end
    vind = Int[]
    vptr = Int[]
    push!(vptr, 1)
    for vertices in vertices_list
        append!(vind, vertices)
        push!(vptr, length(vind)+1)
    end

    #=vind = zeros(Int, length(eind))
    vind[eptr[1:end-1]] .= 1
    vind = cumsum(vind)
    vind = vind[1:end]
    he_prm = sortperm(eind)
    v = eind[he_prm]
    vind = vind[he_prm] 
    vptr = zeros(Int, num_vertices+1)
    println("sizes ", length(vptr), ", ", length(findall(v[1:end-1] .!= v[2:end])))
    vptr[2:end-1] = findall(v[1:end-1] .!= v[2:end]) .+ 1
    vptr[1] = 1
    vptr[end] = eptr[end]=#
    return __hypergraph__(num_vertices, 
                        num_hyperedges, 
                        eptr, 
                        eind, 
                        vptr, 
                        vind, 
                        fixed_vertices, 
                        vwts, 
                        hwts)
end

function remove_single_hyperedges(hgraph::__hypergraph__)
    eptr = hgraph.eptr
    eind = hgraph.eind
    vptr = hgraph.vptr
    vind = hgraph.vind
    new_eptr = [1]
    new_eind = Int[]
    new_hwts = Int[]
    for i in 1:length(eptr)-1
        if eptr[i+1] - eptr[i] > 1
            append!(new_eind, eind[eptr[i]:eptr[i+1]-1])
            push!(new_eptr, length(new_eind)+1)
            push!(new_hwts, hgraph.hwts[i])
        end
    end
    
    return build_hypergraph(hgraph.num_vertices, 
                            length(new_eptr)-1, 
                            new_eptr, 
                            new_eind, 
                            hgraph.fixed, 
                            hgraph.vwts, 
                            new_hwts)
end
