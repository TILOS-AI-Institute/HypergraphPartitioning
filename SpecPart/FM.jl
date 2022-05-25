# FM Related Data Structure Definition

mutable struct FM_Token
    part_vector::Vector{Int}
    cut_size::Int
    part_area::Vector{Int}
    vertices_ordering::Vector{Int}
    vertex_gains::Vector{Int}
    gain_heaps::Vector
    gain_heap_locators::Matrix{Int}
    vertices_part::Vector{Int}
    hyperedges_part::Matrix{Int}
    moved_vertices_flag::Vector{Int}
    marked_vertices_flag::Vector{Int}
    swaps::Vector{Int}
    iter::Int
    pass::Int
end

# Basic Heap Operations

function Heapify!(FM_token::FM_Token, new_vertex_gain::Int, old_vertex_gain::Int, vertex::Int)
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    part_vector = FM_token.part_vector
    part = part_vector[vertex]+1
    heap_location = gain_heap_locators[vertex, 2]
    gain_heaps[part][heap_location][2] = new_vertex_gain

    if old_vertex_gain > new_vertex_gain
        HeapifyDown!(FM_token, new_vertex_gain, heap_location, part)
    else 
        HeapifyUp!(FM_token, heap_location, new_vertex_gain, part)
    end
end

function HeapifyDown!(FM_token::FM_Token, vertex_gain::Int, heap_location::Int, part::Int)
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    vertices_part = FM_token.vertices_part
    max_vertices = vertices_part[part]
    j = heap_location << 1
    gain_left_child = Vector()
    gain_right_child = Vector()
    parent = zeros(Int, 2)
    child = zeros(Int, 2)

    while j <= max_vertices
        gain_left_child = gain_heaps[part][j]

        if j+1 <= max_vertices
            gain_right_child = gain_heaps[part][j+1]

            if gain_left_child[2] < gain_right_child[2]
                j += 1
            end
        else
            if vertex_gain > gain_left_child[2]
                break
            end
        end

        if vertex_gain > gain_left_child[2] && vertex_gain > gain_right_child[2]
            break
        end

        parent = deepcopy(gain_heaps[part][heap_location])
        child = deepcopy(gain_heaps[part][j])
        gain_heaps[part][heap_location] = child
        gain_heap_locators[child[1], 1] = part
        gain_heap_locators[child[1], 2] = heap_location
        gain_heaps[part][j] = parent
        gain_heap_locators[parent[1], 1] = part
        gain_heap_locators[parent[1], 2] = j
        heap_location = j
        j = heap_location << 1
    end
end

function HeapifyUp!(FM_token::FM_Token, heap_location::Int, vertex_gain::Int, part::Int)
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    vertex_gains = FM_token.vertex_gains
    parent = 0

    while heap_location > 1
        parent = heap_location >> 1

        if vertex_gain > gain_heaps[part][parent][2]
            child = gain_heaps[part][heap_location]
            gain_heaps[part][heap_location] = gain_heaps[part][parent]
            heap_node = gain_heaps[part][parent][1]
            gain_heap_locators[heap_node, 1] = part
            gain_heap_locators[heap_node, 2] = heap_location
            gain_heaps[part][parent] = child
            gain_heap_locators[child[1], 1] = part
            gain_heap_locators[child[1], 2] = parent
            heap_location = parent
        else
            break
        end
    end

    #@info "[FM Debug] POST Gain of 25 in heap : $(gain_heaps[gain_heap_locators[25, 1]][gain_heap_locators[25, 2]])"
end

function PopHeap!(FM_token::FM_Token, part::Int)
    success_flag = true
    total_vertices = vertices_part[part]
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    vertices_part = FM_token.vertices_part

    if total_vertices > 1
        gain_block = gain_heaps[part][total_vertices]
        gain_heaps[part][1] = gain_block
        gain_heap_locators[gain_block[1], 1] = part
        gain_heap_locators[gain_block[1], 2] = 1
        HeapifyDown!(FM_token, gain_block[2], 1, part)
        vertices_part[part] -= 1
        success_flag = truevertex_gains
    else
        success_flag = false
    end

    return success_flag
end

function ReturnFromHeap(FM_token::FM_Token, vertex_weights::Vector{Int}, part::Int, max_capacity::Int, min_capacity)
    flip_side_start = false
    heap_location = 1
    best_part = -1
    gain_heaps = FM_token.gain_heaps
    part_area = FM_token.part_area
    vertices_part = FM_token.vertices_part

    if part > -1
        if part_area[part] < min_capacity #|| apart_area[part] > max_capacity
            return (-1)
        end

        best_part = part
    else
        part_0_status = part_area[1] >= min_capacity && part_area[1] <= max_capacity
        part_1_status = part_area[2] >= min_capacity && part_area[2] <= max_capacity

        if part_0_status == false && part_1_status == true
            best_part = 2
        elseif part_0_status == true && part_1_status == false
            best_part = 1
        else
            highest_gain_0 = gain_heaps[1][heap_location][2]
            highest_gain_1 = gain_heaps[2][heap_location][2]
            best_part = highest_gain_0 > highest_gain_1 ? 1 : 2
        end
    end

    best_vertex = gain_heaps[best_part][heap_location][1]
    flip_part = best_part % 2 + 1
    max_depth = 100
    depth = 0

    #@info "[FM Debug] Best vertex before while: $best_vertex from part $best_part"

    while part_area[flip_part] + vertex_weights[best_vertex] > max_capacity && part_area[best_part] - vertex_weights[best_vertex] < min_capacity
        depth += 1
        
        if (depth == max_depth) && (flip_side_start == true)
            return (-1, -1)
        elseif (depth == max_depth) && (flip_side_start == false)
            depth = 1
            best_part = flip_part
            flip_part = best_part % 2 + 1
            heap_location = 1
            flip_side_start = true
        end

        heap_location = heap_location * 2

        if (heap_location > vertices_part[best_part]) && (flip_side_start == false)
            depth = 1
            best_part = flip_part
            flip_part = best_part % 2 + 1
            heap_location = 1
            flip_side_start = true
        elseif heap_location > vertices_part[best_part] && flip_side_start == true
            return (-1, -1)
        end

        best_vertex = gain_heaps[best_part][heap_location][1]
    end

    #@info "[FM Debug] Best vertex after while: $best_vertex from part $best_part"

    return (best_vertex, best_part)
end

# Vertex Move Related Operations

function UpdateAfterMove!(vertex::Int, FM_token::FM_Token, hypergraph_c::Hypergraph_C, incidence_struct::Incidence)
    part_vector = FM_token.part_vector
    hypergraph = hypergraph_c.HG
    fixed_part = hypergraph_c.fixed_part
    i_hedges = incidence_struct.hedges
    i_eptr = incidence_struct.eptr
    h_hedges = hypergraph.hedges
    h_eptr = hypergraph.eptr
    hwts = hypergraph.hwts
    vertices_moved = FM_token.moved_vertices_flag
    vertices_moved[vertex] = FM_token.iter
    vertices_marked = FM_token.marked_vertices_flag
    flip_part = part_vector[vertex]
    part = flip_part == 0 ? 1 : 0
    vertex_gains = FM_token.vertex_gains
    hyperedges_part = FM_token.hyperedges_part
    nbrs = Vector{Int}()

    for i in i_eptr[vertex]:i_eptr[vertex+1]-1
        hyperedge = i_hedges[i]
        hwt = hwts[hyperedge]
        start_idx = h_eptr[hyperedge]
        end_idx = h_eptr[hyperedge+1]
        hyperedges_part[hyperedge, flip_part+1] += 1
        hyperedges_part[hyperedge, part+1] -= 1

        for j in start_idx:end_idx-1
            vertex_j = h_hedges[j]

            if vertex_j == vertex || fixed_part[vertex_j] > -1
                continue
            end

            part_j = part_vector[vertex_j]
            from_part_hyperedge = hyperedges_part[hyperedge, part+1]
            flip_part_hyperedge = hyperedges_part[hyperedge, flip_part+1]

            if part_j == part
                if from_part_hyperedge == 1
                    vertex_gains[vertex_j] += hwt
                end

                if flip_part_hyperedge == 1
                    vertex_gains[vertex_j] += hwt
                end
            else
                if from_part_hyperedge == 0
                    vertex_gains[vertex_j] -= hwt
                end
                
                if flip_part_hyperedge == 2
                    vertex_gains[vertex_j] -= hwt
                end
            end

            if vertices_moved[vertex_j] != FM_token.iter && vertices_marked[vertex_j] != FM_token.pass
                push!(nbrs, vertex_j)
                vertices_marked[vertex_j] = FM_token.pass
            end
        end
    end

    return nbrs
end

function RollBack!(FM_token::FM_Token, hypergraph_c::Hypergraph_C, incidence_struct::Incidence, min_order::Int)
    part_vector = FM_token.part_vector
    part_area = FM_token.part_area
    vertex_weights = hypergraph_c.HG.vwts
    part = -1
    flip_part = -1
    swaps_rev = reverse(FM_token.swaps)

    for i in 1:length(FM_token.swaps)-min_order+1
        vertex = swaps_rev[i]
        part = part_vector[vertex]
        flip_part = part == 0 ? 1 : 0
        part_area[part+1] -= vertex_weights[vertex]
        part_area[flip_part+1] += vertex_weights[vertex]
        part_vector[vertex] = flip_part
        UpdateAfterMove!(vertex, FM_token, hypergraph_c, incidence_struct)
    end
end

function PickBestVertexToMove(FM_token::FM_Token, vertex_weights::Vector{Int}, max_capacity::Int, min_capacity::Int)
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    vertices_part = FM_token.vertices_part    
    heap_location = 0
    best_vertex = -1
    best_part = -1

    if vertices_part[1] == 0 && vertices_part[2] == 0
        return -1
    elseif vertices_part[1] == 0 && vertices_part[2] > 0      
        (best_vertex, best_part) = ReturnFromHeap(FM_token, vertex_weights, 2, max_capacity, min_capacity)
    elseif vertices_part[1] > 0 && vertices_part[2] == 0 
        (best_vertex, best_part) = ReturnFromHeap(FM_token, vertex_weights, 1, max_capacity, min_capacity)
    else
        (best_vertex, best_part) = ReturnFromHeap(FM_token, vertex_weights, -1, max_capacity, min_capacity)
    end

    if best_part > -1
        #@info "[FM Debug] Best Part: $best_part"
        #@info "[FM Debug] Vertices_Part: $(vertices_part)"
        last_gain_block = gain_heaps[best_part][vertices_part[best_part]]
        heap_location = gain_heap_locators[best_vertex, 2]
        gain_heaps[best_part][heap_location] = last_gain_block
        gain_heap_locators[last_gain_block[1], 1] = best_part
        gain_heap_locators[last_gain_block[1], 2] = heap_location
        HeapifyDown!(FM_token, last_gain_block[2], heap_location, best_part)
        vertices_part[best_part] -= 1
    end

    return best_vertex
end

# Heap Manipulation Operations

function UpdateNeighborsInHeap!(FM_token::FM_Token, nbrs::Vector{Int})
    for i in 1:length(nbrs)
        nbr_vtx = nbrs[i]
        new_gain_nbr_vtx = FM_token.vertex_gains[nbr_vtx]
        heap_location_tuple = FM_token.gain_heap_locators[nbr_vtx, :]
        old_gain_nbr_vtx = FM_token.gain_heaps[heap_location_tuple[1]][heap_location_tuple[2]][2]
        Heapify!(FM_token, new_gain_nbr_vtx, old_gain_nbr_vtx, nbr_vtx)
    end
end

function UpdateGainHeap!(FM_token::FM_Token)
    part = 0
    heap_location = 0
    part_vector = FM_token.part_vector
    vertices_ordering = FM_token.vertices_ordering
    vertex_gains = FM_token.vertex_gains
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    vertices_part = FM_token.vertices_part

    for i in 1:length(vertices_ordering)
        vtx = vertices_ordering[i]
        part = part_vector[vtx] + 1
        push!(gain_heaps[part], [vtx, vertex_gains[vtx]])
        vertices_part[part] += 1
        heap_location = vertices_part[part]
        gain_heap_locators[vtx, 1] = part
        gain_heap_locators[vtx, 2] = vertices_part[part]
        HeapifyUp!(FM_token, heap_location, vtx, part)
    end
end

# Data Structure Initialization Operations

function InitializeHyperedges!(part_vector::Vector{Int}, hypergraph::Hypergraph, hyperedges_part::Matrix{Int}, incidence_struct::Incidence; swapped_vertices::Vector{Int} = Int[])
    (start_idx, end_idx) = (0,0)
    hyperedges_list = Vector{Int}()
    hyperedges_flag = zeros(Int, hypergraph.e)

    if isempty(swapped_vertices) == false
        for i in 1:length(swapped_vertices)
            vtx = swapped_vertices[i]
            start_idx = incidence_struct.eptr[vtx]
            end_idx = incidence_struct.eptr[vtx+1]-1
            
            for j in start_idx:end_idx
                hyperedge = incidence_struct.hedges[j]

                if hyperedges_flag[hyperedge] == 0
                    push!(hyperedges_list, hyperedge)
                    hyperedges_flag[hyperedge] = 1
                end
            end
        end
    else
        hyperedges_list = Vector{Int}(1:hypergraph.e)
    end

    for i in 1:length(hyperedges_list)
        hyperedge = hyperedges_list[i]
        start_idx = hypergraph.eptr[hyperedge]
        end_idx = hypergraph.eptr[hyperedge+1]-1

        for j in start_idx:end_idx
            v = hypergraph.hedges[j]
            hyperedges_part[hyperedge, part_vector[v]+1] += 1
        end
    end
end

function InitializeVertices!(part_vector::Vector{Int}, hyperedges_part::Matrix{Int}, h_c::Hypergraph_C, incidence_struct::Incidence, vertex_gains::Vector{Int}; swapped_vertices::Vector{Int} = Int[])
    (part_side, flip_side, start_idx, end_idx, total_vertex_gain, hyperedge_gain, fromDeg, toDeg) = (0,0,0,0,0,0,0,0)
    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    v_hedges = incidence_struct.hedges
    v_eptr = incidence_struct.eptr
    h_eptr = hypergraph.eptr
    hwts = hypergraph.hwts
    vertices_list = Vector{Int}()

    if isempty(swapped_vertices) == false
        vertices_list = swapped_vertices
    else
        vertices_list = Vector{Int}(1:hypergraph.n)
    end

    for i in 1:length(vertices_list)
        vertex = vertices_list[i]
        if fixed_part[vertex] < 0
            part_side = part_vector[vertex] + 1
            flip_side = part_side == 1 ? 2 : 1
            start_idx = v_eptr[vertex]
            end_idx = v_eptr[vertex+1]-1

            for j in start_idx:end_idx
                hyperedge = v_hedges[j]
                hwt = hwts[hyperedge]
                fromDeg = hyperedges_part[hyperedge, part_side]
                toDeg = hyperedges_part[hyperedge, flip_side]

                if toDeg == 0 && fromDeg > 0
                    hyperedge_gain = -hwt
                elseif fromDeg == 1 && toDeg >= 0
                    hyperedge_gain = hwt
                else
                    hyperedge_gain = 0
                end

                total_vertex_gain += hyperedge_gain
            end

            vertex_gains[i] = total_vertex_gain
            total_vertex_gain = 0
        end
    end
end

function ComputeInitialTwoWayGains(part_vector::Vector{Int}, h_c::Hypergraph_C, incidence_struct::Incidence)
    hypergraph = h_c.HG
    e = hypergraph.e
    n = hypergraph.n
    hyperedges_part = zeros(Int, e, 2)
    vertex_gains = zeros(Int, n)
    InitializeHyperedges!(part_vector, hypergraph, hyperedges_part, incidence_struct)
    InitializeVertices!(part_vector, hyperedges_part, h_c, incidence_struct, vertex_gains)

    return (hyperedges_part, vertex_gains)
end

function InitializeFM(h_c::Hypergraph_C, incidence_struct::Incidence, part_vector::Vector{Int}, part_area::Vector{Int})
    hypergraph = h_c.HG
    n = hypergraph.n
    gain_heaps = [Vector() for i in 1:2]
    gain_heap_locators = zeros(Int, n, 2)
    moved_vertices_flag = zeros(Int, n)
    marked_vertices_flag = zeros(Int, n)
    vertices_part = zeros(Int, 2)
    (cut_size, part_area) = FindCutSize(part_vector, hypergraph) 
    (hyperedges_part, vertex_gains) = ComputeInitialTwoWayGains(part_vector, h_c, incidence_struct)
    vertices_ordering = sortperm(vertex_gains, rev=true)

    return FM_Token(part_vector, cut_size, part_area, vertices_ordering, vertex_gains, gain_heaps, gain_heap_locators, vertices_part, hyperedges_part, moved_vertices_flag, marked_vertices_flag, Int[], 0, 0)
end

function ReInitializeFM(FM_token::FM_Token)
    part = 0
    heap_location = 0
    vertices_part = FM_token.vertices_part
    vertices_part[1] = 0
    vertices_part[2] = 0
    part_vector = FM_token.part_vector
    vertex_gains = FM_token.vertex_gains
    gain_heaps = FM_token.gain_heaps
    gain_heap_locators = FM_token.gain_heap_locators
    vertices_ordering = sortperm(vertex_gains, rev=true)
    FM_token.vertices_ordering = vertices_ordering
    FM_token.swaps = Int[]
    FM_token.iter = 0
    FM_token.pass = 0

    for i in 1:length(vertices_ordering)
        vtx = vertices_ordering[i]
        part = part_vector[vtx] + 1
        vertices_part[part] += 1
        
        if vertices_part[part] > length(gain_heaps[part])
            push!(gain_heaps[part], [vtx, vertex_gains[vtx]])
        else
            gain_heaps[part][vertices_part[part]] = [vtx, vertex_gains[vtx]]
        end

        heap_location = vertices_part[part]
        gain_heap_locators[vtx, 1] = part
        gain_heap_locators[vtx, 2] = vertices_part[part]
        HeapifyUp!(FM_token, heap_location, vtx, part)
    end
end

function FindCutSize(part_vector::Vector{Int}, hypergraph::Hypergraph)
    hwts = hypergraph.hwts
    vwts = hypergraph.vwts
    cut_size = 0
    part_area = zeros(Int, 2)

    for i in 1:hypergraph.e
        start_idx = hypergraph.eptr[i]
        end_idx = hypergraph.eptr[i+1]-1
        part_hash = 0

        for j in start_idx:end_idx
            v = hypergraph.hedges[j]
            part_hash += part_vector[v]
        end

        if part_hash > 0 && part_hash < (end_idx-start_idx + 1)
            cut_size += hwts[i]
        end
    end

    for i in 1:hypergraph.n
        part_area[part_vector[i]+1] += vwts[i]
    end

    return (cut_size, part_area)
end

# FM Driver Operation

function FM(h_c::Hypergraph_C, incidence_struct::Incidence, part_vector::Vector{Int}, part_area::Vector{Int}, max_capacity::Int, min_capacity::Int; fm_max_iters::Int = 10, fm_max_passes::Int = 100, bad_moves_threshold::Int = 100, early_exit::Bool = true)
    if early_exit == false
        fm_max_passes = h_c.HG.n
        bad_moves_threshold = h_c.HG.n
    end
    
    saturation_depth = 0
    max_saturation_depth = 2
    min_pos = 0
    bad_moves = 0
    vertex_weights = h_c.HG.vwts
    FM_token = InitializeFM(h_c, incidence_struct, part_vector, part_area)
    cut_chain = Vector{Int}([FM_token.cut_size])
    local_min_cut_size = FM_token.cut_size
    global_min_cut_size = local_min_cut_size
    prev_local_min_cut_size = local_min_cut_size
    global_part_area = FM_token.part_area
    global_part_vector = zeros(Int, h_c.HG.n)
    UpdateGainHeap!(FM_token)

    for i in 1:fm_max_iters
        if saturation_depth == max_saturation_depth
            break
        end

        min_pos = 0
        FM_token.swaps = Vector{Int}()
        FM_token.iter = i

        for j in 1:fm_max_passes
            FM_token.pass = j
            best_vertex = PickBestVertexToMove(FM_token, vertex_weights, max_capacity, min_capacity)
            
            if best_vertex == -1
                break
            end

            push!(FM_token.swaps, best_vertex)
            from = FM_token.part_vector[best_vertex]
            to = from == 0 ? 1 : 0
            FM_token.part_vector[best_vertex] = to
            FM_token.part_area[from+1] -= vertex_weights[best_vertex]
            FM_token.part_area[to+1] += vertex_weights[best_vertex]
            FM_token.cut_size -= FM_token.vertex_gains[best_vertex]

            if FM_token.cut_size < local_min_cut_size
                local_min_cut_size = FM_token.cut_size
                min_pos = j
            end

            #@info "[FM Debug] Iter $i :: Pass $j :: Best vertex $best_vertex :: Gain $(FM_token.vertex_gains[best_vertex]) :: cut $(FM_token.cut_size) :: area $(FM_token.part_area)"

            push!(cut_chain, FM_token.cut_size)
            neighbors_best_vertex = UpdateAfterMove!(best_vertex, FM_token, h_c, incidence_struct)
            UpdateNeighborsInHeap!(FM_token, neighbors_best_vertex)
            bad_moves = cut_chain[end] - cut_chain[end-1] > 0 ? bad_moves + 1 : 0

            if bad_moves >= bad_moves_threshold 
                break
            end
        end

        #@info "[FM Debug] FM Swaps: $(FM_token.swaps)"
        #@info "[FM Debug] Min Pos: $min_pos"

        if length(FM_token.swaps) > 0
            FM_token.cut_size = local_min_cut_size
            RollBack!(FM_token, h_c, incidence_struct, min_pos+1)
            FM_token.swaps = min_pos > 0 ? FM_token.swaps[1:min_pos] : Int[]
            rollback_vertex_gains = min_pos > 0 ? zeros(Int, length(FM_token.swaps)) : zeros(Int, h_c.HG.n)

            @info "[FM] ITERATION $i :: CUT RECORDED: $(FM_token.cut_size) :: AREA : $(FM_token.part_area)"

            InitializeVertices!(FM_token.part_vector, FM_token.hyperedges_part, h_c, incidence_struct, rollback_vertex_gains, swapped_vertices = FM_token.swaps)
        
            if min_pos > 0
                for i in 1:length(rollback_vertex_gains)
                    vertex = FM_token.swaps[i]
                    vertex_gain = rollback_vertex_gains[i]
                    FM_token.vertex_gains[vertex] = vertex_gain
                end
            end

            ReInitializeFM(FM_token)
        else
            @info  "[FM] ITERATION $i FAILED TO FIND FEASIBLE VERTEX TO MOVE!! EXITING FM!"
            break
        end

        saturation_depth = local_min_cut_size >= prev_local_min_cut_size ? saturation_depth + 1 : 0

        if local_min_cut_size < global_min_cut_size
            global_min_cut_size = local_min_cut_size
            global_part_vector = FM_token.part_vector
            global_part_area = FM_token.part_area
        end
    end

    return (global_part_vector, global_min_cut_size, global_part_area)
end