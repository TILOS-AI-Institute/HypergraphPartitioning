struct __hypergraph__
    num_vertices::Int
    num_hyperedges::Int
    eptr::Vector{Int}
    eind::Vector{Int}
    vptr::Vector{Int}
    vind::Vector{Int}
    fixed::Vector{Int}
    vwts::Vector{Int}
    hwts::Vector{Int}
end

struct __pindex__
    p1::Vector{Int}
    p2::Vector{Int}
end

struct __least_common_ancestor__
    rmq_sparse_table::Matrix{Int}
    euler_level::Vector{Int}
    child::Vector{Int}
    parents::Vector{Int}
    euler_tour::Vector{Int}
    level_vec::Vector{Int}
    fts::Vector{Int}
    ifts::Vector{Int}
end

struct __cut_profile__
    vtx_cuts::AbstractArray{Int}
    edge_cuts::Vector{Int}
    edge_diff::Vector{Int}
    pred::Vector{Int}
    edge_terminators::Vector{Int}
    p::__pindex__
    forced_type::Vector{Int}
    forced_0::Vector{Int}
    forced_1::Vector{Int}
    forced_01::Vector{Int}
    FB0::Vector{Int}
    FB1::Vector{Int}
    edge_cuts_0::Vector{Int}
    edge_cuts_1::Vector{Int}
end

struct __best_partition__
    total_cost::Float64
    area_cost::Float64
    cut_cost::Float64
    cutsize::Int
    partition::Vector{Int}
    cut_point::Int
    area::Vector{Int}
end

struct __recursive_parts__
    hypergraph::__hypergraph__
    T::SimpleWeightedGraphs.SimpleGraph
    distilled_cuts::__cut_profile__
    capacities::Vector{Int}
    cluster_labels::Vector{Int}
end

mutable struct __tree_cuts__
    nforced0::Int
    nforced1::Int
    nforced01::Int
    total_vwts::Int
    exc0::Vector{Int}
    exc1::Vector{Int}
    area::Vector{Int}
    cut_cost0::Vector{Float64}
    cut_cost1::Vector{Float64}
    ratio_cost::Vector{Float64}
    cut_cost::Vector{Float64}
    area_cost::Vector{Float64}
    total_cost::Vector{Float64}
    polarity::Vector{Int}
    status_flag::Vector{Int}
    area_util0::Vector{Int}
    area_util1::Vector{Int}
    pred::Vector{Int}
    hyperedges_flag::Vector{Int}
end

function hypl(hypergraph::__hypergraph__, x::AbstractArray, epsilon::Int)
    eind = hypergraph.eind
    eptr = hypergraph.eptr
    n = length(x)
    m = hypergraph.num_hyperedges
    y = zeros(Float64, n)
    w = hypergraph.hwts

    for j in 1:m
        first_valid_entry = eptr[j]
        first_invalid_entry = eptr[j+1]
        k = first_invalid_entry - first_valid_entry
        scale = (floor(k/2) * ceil(k/2))/(k-1)
        sm = 0.0
        for t in first_valid_entry:first_invalid_entry-1
            sm += x[eind[t]]
        end
        sm /= k
        for t in first_valid_entry:first_invalid_entry-1
            idx = eind[t]
            #y[idx] += w[j] * (x[idx] - sm)/scale
            y[idx] += w[j] * (x[idx] - sm)/(scale*epsilon)
        end
    end
    return y
end

function clique(x::AbstractArray, 
                vwts::Vector{Int}, 
                multiplier::AbstractArray)
    twt = sum(vwts)
    n = size(x, 1)
    y = zeros(Float64, n)
    s = multiplier[1]/twt
    kvec = vwts'x
    @sync Threads.@threads for j in 1:n
        y[j] += twt * ((vwts[j] * x[j]) - ((kvec * vwts[j])/twt)) * s
    end
    return y
end

function bi_clique(x::AbstractArray, pindex::__pindex__)
    n = length(x)
    y = zeros(n)
    d1 = ones(length(pindex.p1))
    d2 = ones(length(pindex.p2))
    t1 = Threads.@spawn (sum(d2) .* d1 .* x[pindex.p1] - d1 * (d2' * x[pindex.p2]))
    t2 = Threads.@spawn (sum(d1) .* d2 .* x[pindex.p2] - d2 * (d1' * x[pindex.p1])) 
    t1 = fetch(t1)
    t2 = fetch(t2)
    y[pindex.p1] = t1
    y[pindex.p2] = t2
    return y
end

function process_hint(partition::Vector{Int}, 
                    new_indices::Vector{Int},
                    processed_partition::Vector{Int})
    for i in 1:length(new_indices)
        mapped_i = new_indices[i]
        if mapped_i == 0
            continue
        end
        processed_partition[mapped_i] = partition[i]
    end
end 

function find_labels(clusters::AbstractArray, 
                    n::Int)
    labels = zeros(Int, n)
    for i in 1:length(clusters) 
        for j in 1:length(clusters[i])
            labels[clusters[i][j]] = i-1
        end
    end
    return labels
end

function write_partition(partition::Vector{Int},
                        partition_file_name::String)
    f = open(partition_file_name, "w")
    for i in 1:length(partition) 
        println(f, partition[i])
    end
    close(f)
end

metis_path = "/home/bodhi91/SpecPart_extension/metis"