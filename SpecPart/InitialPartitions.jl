function GenerateRandomInitialPartitions(vertex_weights::Vector{Int}, max_capacity::Int, min_capacity::Int)
    n = length(vertex_weights)
    vertices_list = Vector{Int}(1:n)
    Shuffle.shuffle!(vertices_list)
    part_vector = -ones(Int, n)
    part_area = zeros(Int, 2)

    for i in 1:length(vertices_list)
        v = vertices_list[i]
        possible_area = part_area[1] + vertex_weights[v] 

        if possible_area < max_capacity
            part_vector[v] = 0
            part_area[1] = possible_area
        else
            part_vector[v] = 1
            part_area[2] += vertex_weights[v]
        end
    end

    return (part_vector, part_area)
end
