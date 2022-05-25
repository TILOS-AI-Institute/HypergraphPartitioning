mutable struct Hypergraph
	n::Int
	e::Int
	hedges::Vector{Int}
	eptr::Vector{Int}
    vwts::Vector{Int}
    hwts::Vector{Int}
end

mutable struct Hypergraph_C
    HG::Hypergraph
    incidence_list::Vector{Vector{Int}}
    hyperedges_pair_list::Array{Int}
    community::Vector{Int}
    fixed_part::Vector{Int}
    fixed_vertex_flag::Bool
    hyperedges_hash::Dict{Int, Int}
end

struct Incidence
	n::Int
	e::Int
	hedges::Vector{Int}
	eptr::Vector{Int}
	d::Vector{Int}
end
