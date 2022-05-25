using LightGraphs
using SimpleTraits
using SimpleWeightedGraphs
using LDLFactorizations
using Statistics
using IterativeSolvers
using Laplacians
using LinearAlgebra
using LinearMaps
using DataFrames
using CSV
using Clustering
using ParallelKMeans
using Base: Float64, show_enclosed_list

struct LCA
    rmqSparseTable::Matrix{Int}
    eulerLevel::Vector{Int}
    child::Vector{Int}
    parents::Vector{Int}
    eulerTour::Vector{Int}
    levelVec::Vector{Int}
    fts::Vector{Int}
    ifts::Vector{Int}
end

struct Pindex
    p1::Vector{Int}
    p2::Vector{Int}
end

struct Instance
    H::Hypergraph
    A::SparseMatrixCSC
    B::Incidence
    p::Pindex
end

struct CutProfile
    vtxCuts::AbstractArray{Int}
    edgeCuts::Vector{Int}
    edgeDiff::Vector{Int}
    pred::Vector{Int}
    edgeTerminators::Vector{Int}
    pindex::Pindex
    forcedType::Vector{Int}
    forced_0::Vector{Int}
    forced_1::Vector{Int}
    forced_01::Vector{Int}
    FB0::Vector{Int}
    FB1::Vector{Int}
    edgeCuts0::Vector{Int}
    edgeCuts1::Vector{Int}
end

struct BestCut
    tCost::Float64
    areaCost::Float64
    cutCost::Float64
    cut_size::Int64 
    bins::Vector{Int64}
    cutpoint::Int64
    areaPart::Vector{Int}
end 

include("GenTrees.jl")
include("rmq/EulerTour.jl")
include("rmq/lca2rmq.jl")
include("rmq/NodeLevels.jl")
include("rmq/Queries.jl")
include("rmq/Queries.jl")
include("rmq/rmq.jl")
include("rmq/rmq_solve.jl")
#include("SweepCuts.jl")
include("TreeSweepUtilities.jl")
include("TreeBipartition.jl")
include("cmg/src/cmgAlg.jl")
include("cmg/src/CombinatorialMultigrid.jl")
include("hypL.jl")
include("dCliques.jl")
include("dBiClique.jl")
include("SolveEigenVecs.jl")
include("HypergraphToGraph.jl")
include("CutProfile.jl")
include("Clustering.jl")
include("IsolateIslands.jl")
include("Hyperstretch.jl")
include("GenerateEigenAwareGraph.jl")
include("PlotEmbedding.jl")
include("Modified_Prims.jl")
include("BuildMetisGraph.jl")
include("AnalyzeSpectralCuts.jl")
include("IterativeMSTReconstruction.jl")

function metis(fname::String)
    log_name = "metis" * string(rand(1:1000)) * "log.txt"
    #metis_r = "./gpmetis " * fname * " 2 -ptype=rb -ufactor=100 -dbglvl=0" # > 'metis_log.txt'"""
    #metis_r = `./gpmetis $fname 2 -ptype=rb -ufactor=100 -dbglvl=0` 
    metis_r = `./metis_script.sh $fname $log_name`
    run(metis_r, wait=true)
    cmd = "rm -r " * log_name
    run(`sh -c $cmd`, wait=true)
end

function WriteTreeToFile(tree::SimpleWeightedGraph, fname::String)
    f = open(fname, "w")
    
    println(f, nv(tree), " ", ne(tree))

    for e in edges(tree)
        println(f, e.src, ",", e.dst, " ", e.weight)
    end
    
    close(f)
end

function ComputeTreePartition(tree_eig::SimpleWeightedGraph, X::Array{Float64}, wt_matrix::SparseMatrixCSC, hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vertices::Pindex, tree_type::Int)
    tree, treeMatrix = GenEigenTree(tree_eig, wt_matrix, X, tree_type)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    return pvec
end

function GenerateParallelCliques(adj_mat::SparseMatrixCSC, X::Array{Float64})
    clique_tasks = Vector{Task}()
    clique_graphs = [SimpleWeightedGraph() for i in 1:6]
    push!(clique_tasks, Threads.@spawn ModifyExpanderWts(adj_mat, X[:, 1], 1, false))
    push!(clique_tasks, Threads.@spawn ModifyExpanderWts(adj_mat, X[:, 2], 1, false))
    push!(clique_tasks, Threads.@spawn ModifyExpanderWts(adj_mat, X, 2, false))
    push!(clique_tasks, Threads.@spawn ModifyExpanderWts(adj_mat, X[:, 1], 1, true))
    push!(clique_tasks, Threads.@spawn ModifyExpanderWts(adj_mat, X[:, 2], 1, true))
    push!(clique_tasks, Threads.@spawn ModifyExpanderWts(adj_mat, X, 2, true))

    for i in 1:length(clique_tasks)
        clique_graphs[i] = fetch(clique_tasks[i])
    end

    return clique_graphs
end

function GenerateParallelPartitions(clique_graphs::Vector{SimpleWeightedGraph{Int, Float64}}, X::Array{Float64}, wt_matrix::SparseMatrixCSC, hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vertices::Pindex)
    partitions_tasks = Vector{Task}()
    partition_tokens = [Vector{Int}() for i in 1:8]

    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[1], X[:, 1], wt_matrix, hypergraph, incidence_struct, fixed_vertices, 2))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[2], X[:, 2], wt_matrix, hypergraph, incidence_struct, fixed_vertices, 2))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[3], X, wt_matrix, hypergraph, incidence_struct, fixed_vertices, 2))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[1], X[:, 1], wt_matrix, hypergraph, incidence_struct, fixed_vertices, 4))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[2], X[:, 2], wt_matrix, hypergraph, incidence_struct, fixed_vertices, 4))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[4], X[:, 1], wt_matrix, hypergraph, incidence_struct, fixed_vertices, 1))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[5], X[:, 2], wt_matrix, hypergraph, incidence_struct, fixed_vertices, 1))
    push!(partitions_tasks, Threads.@spawn ComputeTreePartition(clique_graphs[6], X, wt_matrix, hypergraph, incidence_struct, fixed_vertices, 1))
    
    for i in 1:length(partitions_tasks)
        partition_tokens[i] = fetch(partitions_tasks[i])
    end

    return partition_tokens
end

function GenerateFamilyOfTrees(h_c::Hypergraph_C, incidence_struct::Incidence, ubfactor::Int; eigenvecs::Int = 1, expander_cycles::Int = 3, eigen_iters::Int = 50, bestsolns_on_tree::Int = 3, families::Int = 1)
    @info "FINDING SPECTRAL CUT"
    line_log = repeat("=", 80)
    @info "$line_log"

    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    max_capacity = Int(ceil(sum(hypergraph.vwts) * (50+ubfactor)/100))
    capacities = [max_capacity, max_capacity]

    binMatrix = zeros(Int, families*8, hypergraph.n)
    list_of_cuts = zeros(Int, families*8)
    cutsizes = zeros(Int, families, 8)
    pcuts = zeros(Int, 8, hypergraph.n)
    ptr = 0

    expander_cycles = [1] #, 3, 5]

    #GenerateParallelAdjMats()


    (adj_mat, wt_matrix) = HypergraphToGraph(hypergraph, 1)
    X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)

    pvec_threads = []
    pvecs = Vector{Vector{Int}}()

    for i in 1:families
        clique_graphs = GenerateParallelCliques(adj_mat, X)
        pvecs = GenerateParallelPartitions(clique_graphs, X, wt_matrix, hypergraph, incidence_struct, fixed_vertices)
        
        for i in 1:length(pvecs)
            ptr += 1
            binMatrix[ptr, :] = pvecs[i]
            cut_size, ~ = FindCutSize(pvecs[i], hypergraph)
            @info "[ITER $ptr] :: CUT RECORDED :: $cut_size"
            list_of_cuts[ptr] = cut_size
        end

        #(adj_mat, wt_matrix) = HypergraphToGraph(hypergraph, expander_cycles[i])
    end

    p_perm = sortperm(list_of_cuts)
    

    #=for i in 1:families
        ptr += 1
        tree_eig_1 = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
        tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 2)
        pcuts[1, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[1, :], cutsizes[i, 1] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[1, :]

        ptr += 1
        tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
        tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 2)
        pcuts[2, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[2, :], cutsizes[i, 2] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[2, :]

        ptr += 1
        tree_eig_12 = ModifyExpanderWts(adj_mat, X, eigenvecs, false)
        tree, treeMatrix = GenEigenTree(tree_eig_12, wt_matrix, X, 2)
        pcuts[3, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[3, :], cutsizes[i, 3] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[3, :]

        ptr += 1
        tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 4)
        pcuts[4, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[4, :], cutsizes[i, 4] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[4, :]

        ptr += 1
        tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 4)
        pcuts[5, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[5, :], cutsizes[i, 5] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[5, :]

        ptr += 1
        tree_eig_1 = ModifyExpanderWts(adj_mat, X[:, 1], 1, true)
        tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 1)
        pcuts[6, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[6, :], cutsizes[i, 6] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[6, :]

        ptr += 1
        tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, true)
        tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 1)
        pcuts[7, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[7, :], cutsizes[i, 7] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[7, :]

        ptr += 1
        tree_eig_12 = ModifyExpanderWts(adj_mat, X, eigenvecs, true)
        tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X, 1)
        pcuts[8, :] = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
        #pcuts[8, :], cutsizes[i, 8] = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        binMatrix[ptr, :] = pcuts[8, :]

        #mincut, mintree = findmin(cutsizes[i, :])
        #@info "[ITER $i] :: MIN CUT RECORDED :: $mincut"
        #binMatrix[i, :] = pcuts[mintree, :] .- 1
        
        (adj_mat, wt_matrix) = HypergraphToGraph(hypergraph, expander_cycles[i])
    end=#

    return binMatrix[p_perm[1:8], :]
end

function GenerateCutFromTreeUsingSpectral(tree::SimpleGraph, tree_matrix::SparseMatrixCSC, wt_matrix::SparseMatrixCSC, h_c::Hypergraph_C, incidence_struct::Incidence, capacities::Vector{Int}, eigenvecs::Int, eigen_iters::Int)
    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    
    num_vertices = nv(tree)
    num_hyperedges = ne(tree)
    hedges = Vector{Int}()
    eptr = zeros(Int, num_hyperedges+1)
    hwts = ones(Int, num_hyperedges)
    eptr[1] = 1
    l = 1
    i = 0

    for e in edges(tree)
        i += 1
        push!(hedges, e.src)
        push!(hedges, e.dst)
        l += 2
        eptr[i] = l
        hwts[i] = wt_matrix[e.src, e.dst]
    end

    hypergraph_tree = Hypergraph(num_vertices, num_hyperedges, hedges, eptr, hypergraph.vwts, hwts)

    ExportHypergraph(hypergraph_tree, "tree_hypergraph.hgr", -ones(Int, num_vertices))

    hmetis("tree_hypergraph.hgr", 2, 5, 10, 1, 1, 0, 1, 24)

    f = open("tree_hypergraph.hgr.part.2")

    bin = zeros(Int, num_vertices)
    i = 0

    for ln in eachline(f)
        i += 1
        bin[i] = parse(Int, ln)
    end

    cut, area = FindCutSize(bin, hypergraph)

    println("Cut from hmetis: ", cut, ", area split: ", area)

    return bin

    X = GenEigenVecs(hypergraph_tree, tree_matrix, fixed_vertices, false, eigenvecs, eigen_iters)
    binMatrix = BestCutsOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 3)
    return binMatrix
end

function METIS_Run(hypergraph::Hypergraph, tree::SimpleWeightedGraph)
    fname = BuildMetisGraph(tree)
    metis(fname)

    f = open(fname*".part.2", "r")
    i = 0
    pvec = zeros(Int, hypergraph.n)

    for ln in eachline(f)
        i += 1
        pvec[i] = parse(Int, ln)
    end
    
    close(f)

    cut_size = FindCutSize(pvec, hypergraph), pvec
    cmd = "rm -r " * fname
    run(`sh -c $cmd`, wait=true)
    cmd = "rm -r " * fname * ".part.2"
    run(`sh -c $cmd`, wait=true)
    
    
    return cut_size
end

function MetisOnTree(hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vertices::Pindex, tree::SimpleGraph, treeMatrix::SparseMatrixCSC)
    cutInfo = analyzeCutsOnTree(hypergraph, incidence_struct, hypergraph.vwts, fixed_vertices, tree, ones(Int, hypergraph.e))
    forced_0 = cutInfo.forced_0
	forced_1 = cutInfo.forced_1
	forced_01 = cutInfo.forced_01
    FB0 = cutInfo.FB0
	FB1 = cutInfo.FB1
	edgeCuts0 = cutInfo.edgeCuts0
	edgeCuts1 = cutInfo.edgeCuts1
    edgeCuts = cutInfo.edgeCuts
    nforced_0 = sum(hypergraph.hwts[forced_0])
    nforced_1 = sum(hypergraph.hwts[forced_1])
    nforced_01 = sum(hypergraph.hwts[forced_01])
    exc0 = zeros(Int, hypergraph.n)
    exc1 = zeros(Int, hypergraph.n)
    cutCost = similar(exc0)

    for i in 1:hypergraph.n
        if cutInfo.pred[i] == i
            cutCost[i] = 1e09
            continue
        end
            
        exc0[i] = edgeCuts[i] + nforced_0 - FB0[i] + edgeCuts1[i] + nforced_01
		exc1[i] = edgeCuts[i] + nforced_1 - FB1[i] + edgeCuts0[i] + nforced_01
        cutCost[i] = min(exc0[i], exc1[i])
        i_parent = cutInfo.pred[i]
        treeMatrix[i, i_parent] = cutCost[i]
        treeMatrix[i_parent, i] = cutCost[i]
    end

    g = SimpleWeightedGraph(treeMatrix)
    cut_elements, pvec = METIS_Run(hypergraph, g)

    #@info "METIS CUT FOUND IS: $(cut_elements[1]) WITH SPLIT: $(cut_elements[2])"

    return pvec #, cut_elements[1]
end

function PartTreeUsingMetis(h_c::Hypergraph_C, incidence_struct::Incidence, ubfactor::Int; eigenvecs::Int = 1, expander_cycles::Int = 1, eigen_iters::Int = 50, bestsolns_on_tree::Int = 3)
    @info "FINDING SPECTRAL CUT"
    line_log = repeat("=", 80)
    @info "$line_log"

    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    max_capacity = Int(ceil(sum(hypergraph.vwts) * (50+ubfactor)/100))
    capacities = [max_capacity, max_capacity]

    binMatrix = zeros(Int, 8, hypergraph.n)

    (adj_mat, wt_matrix) = HypergraphToGraph(hypergraph, expander_cycles)
    X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)

    tree_eig_1 = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
    tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 2)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[1, :] = pvec

    tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 2)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[2, :] = pvec

    tree_eig_12 = ModifyExpanderWts(adj_mat, X, eigenvecs, false)
    tree, treeMatrix = GenEigenTree(tree_eig_12, wt_matrix, X, 2)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[3, :] = pvec

    tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 4)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[4, :] = pvec

    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 4)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[5, :] = pvec

    tree_eig_1 = ModifyExpanderWts(adj_mat, X[:, 1], 1, true)
    tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 1)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[6, :] = pvec

    tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, true)
    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 1)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[7, :] = pvec

    tree_eig_12 = ModifyExpanderWts(adj_mat, X, 2, true)
    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X, 1)
    pvec = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix)
    binMatrix[8, :] = pvec

    return binMatrix
end

function GenerateSpectralCut(h_c::Hypergraph_C, incidence_struct::Incidence, ubfactor::Int; eigenvecs::Int = 1, expander_cycles::Int = 1, eigen_iters::Int = 50, bestsolns_on_tree::Int = 3)
    @info "FINDING SPECTRAL CUT"
    line_log = repeat("=", 80)
    @info "$line_log"

    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    max_capacity = Int(ceil(sum(hypergraph.vwts) * (50+ubfactor)/100))
    capacities = [max_capacity, max_capacity]

    (adj_mat, wt_matrix) = HypergraphToGraph(hypergraph, expander_cycles)

    #adj_mat_graph = SimpleWeightedGraph(adj_mat)
    #WriteTreeToFile(adj_mat_graph, "expander_graph.dat")

    X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)

    #=tree_eig = ModifyExpanderWts(adj_mat, X, eigenvecs, false)
    tree, treeMatrix = GenEigenTree(tree_eig, wt_matrix, X, 2)
    
    for i in 1:3    
        tree_cut_token = TreeCutToken(0, 0, 0, 0, zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), zeros(Int, 2), zeros(Float64, hypergraph.n),  zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n),
                        zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n), zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), 
                        zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), ones(Int, hypergraph.e))
        feasible_cuts, feasible_edges = FindFeasibleEdges(tree, tree_cut_token, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        if isempty(feasible_cuts)
            @info "[ITER $i] NO FEASIBLE CUT FOUND!"
        else
            @info "FEASIBLE CUTS: $feasible_cuts"
            @info "[ITER $i] BEST CUT: $(minimum(feasible_cuts))"
        end


        IterativeMSTReconstruction(tree_eig, feasible_edges)
        tree, ~ = GenEigenTree(tree_eig, wt_matrix, X, 2)
    end

    return X=#

    tree_family = []

    # Using MST 

    tree_eig_1 = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
    tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 2)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"
    push!(tree_family, tree)

    tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 2)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"
    push!(tree_family, tree)

    tree_eig_12 = ModifyExpanderWts(adj_mat, X, eigenvecs, false)
    tree, treeMatrix = GenEigenTree(tree_eig_12, wt_matrix, X, 2)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"
    push!(tree_family, tree)

    # Using Path

    tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 4)
    push!(tree_family, tree)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"

    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 4)
    push!(tree_family, tree)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"

    # Using LST

    tree_eig_1 = ModifyExpanderWts(adj_mat, X[:, 1], 1, true)
    tree, treeMatrix = GenEigenTree(tree_eig_1, wt_matrix, X[:, 1], 1)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"
    push!(tree_family, tree)

    tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, true)
    tree, treeMatrix = GenEigenTree(tree_eig_2, wt_matrix, X[:, 2], 1)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"
    push!(tree_family, tree)

    tree_eig_12 = ModifyExpanderWts(adj_mat, X, eigenvecs, true)
    tree, treeMatrix = GenEigenTree(tree_eig_12, wt_matrix, X, 1)
    #cutsize, part_area = METIS_Run(hypergraph, SimpleWeightedGraph(treeMatrix))
    #@info "METIS CUT SIZE: $cutsize WITH SPLIT $part_area"
    push!(tree_family, tree)

    #=path1 = sortperm(X[:,1])
    path2 = sortperm(X[:,2])
    roots = Int[]
    middle = Int(round(hypergraph.n/2))

    append!(roots, path1[middle-2:middle+2])
    append!(roots, path2[middle-2:middle+2])

    tree_eig_2 = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
    short_trees = GenerateShortestPathTrees(tree_eig_2, roots)
    
    append!(tree_family, short_trees)=#

    bin_family = zeros(Int, length(tree_family), hypergraph.n)
    min_cut = 1e12
    best_tree = -1
    dim = -1

    for i in 1:length(tree_family)
        if i < 4
            @info "MST Cuts: \n"
        elseif i < 6
            @info "PATH Cut: \n"
        elseif i < 9
            @info "LST Cuts: \n"
        else
            @info "Shorted Path Cuts: \n"
        end

        bin_family[i, :], cut_size = FindBestCutOnTree(tree_family[i], hypergraph, incidence_struct, fixed_vertices, capacities, 1)
        #push!(bin_family, BestCutsOnTree(tree_family[i], hypergraph, incidence_struct, fixed_vertices, capacities, 1))

        if cut_size < min_cut
            min_cut = cut_size
            best_tree = i
            if i == 1 || i == 6
                dim = 1
            elseif i == 2 || i == 7
                dim = 2
            elseif i == 3 || i == 8
                dim = 0
            end
        end
    end

    tree_type = 0

    if best_tree <= 3
        tree_type = 2
    elseif 3 < best_tree <= 5
        tree_type = 4
    else
        tree_type = 1
    end


    return bin_family[best_tree,: ], dim, tree_type

    #=

    # Using MST

    tree_eig = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
    tree1, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 1], 2)

    tree_eig = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
    tree2, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 2], 2)

    tree_eig = ModifyExpanderWts(adj_mat, X, eigenvecs, false)
    tree3, ~ = GenEigenTree(tree_eig, wt_matrix, X, 2)

    # Using Modified MST

    tree_eig = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
    tree8, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 1], 3)

    tree_eig = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
    tree9, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 2], 3)

    tree_eig = ModifyExpanderWts(adj_mat, X, eigenvecs, false)
    tree10, ~ = GenEigenTree(tree_eig, wt_matrix, X, 3)

    # Using Path

    tree4, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 1], 4)

    # Using LST

    tree_eig = ModifyExpanderWts(adj_mat, X[:, 1], 1, true)
    tree5, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 1], 1)

    tree_eig = ModifyExpanderWts(adj_mat, X[:, 2], 1, true)
    tree6, ~ = GenEigenTree(tree_eig, wt_matrix, X[:, 2], 1)

    tree_eig = ModifyExpanderWts(adj_mat, X, eigenvecs, true)
    tree7, ~ = GenEigenTree(tree_eig, wt_matrix, X, 1)

    #roots = []
    #trees = GenerateShortestPathTrees(tree_eig, roots)
    

    @info "MST Cuts: \n"
    bin1 = BestCutsOnTree(tree1, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    bin2 = BestCutsOnTree(tree2, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    bin3 = BestCutsOnTree(tree3, hypergraph, incidence_struct, fixed_vertices, capacities, 1)

    #@info "MMST Cuts: \n"
    #bin8 = BestCutsOnTree(tree8, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    #bin9 = BestCutsOnTree(tree9, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    #bin10 = BestCutsOnTree(tree10, hypergraph, incidence_struct, fixed_vertices, capacities, 1)

    @info "Path Cut: \n"
    bin4 = BestCutsOnTree(tree4, hypergraph, incidence_struct, fixed_vertices, capacities, 1)

    @info "LST Cuts: \n"
    bin5 = BestCutsOnTree(tree5, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    bin6 = BestCutsOnTree(tree6, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    bin7 = BestCutsOnTree(tree7, hypergraph, incidence_struct, fixed_vertices, capacities, 1)
    
    binMatrix = zeros(Int, 7, hypergraph.n)
    binMatrix[1, :] = bin1[1, :]
    binMatrix[2, :] = bin2[1, :]
    binMatrix[3, :] = bin3[1, :]
    binMatrix[4, :] = bin4[1, :]
    binMatrix[5, :] = bin5[1, :]
    binMatrix[6, :] = bin6[1, :]
    binMatrix[7, :] = bin7[1, :]
    #binMatrix[8, :] = bin8[1, :]
    #binMatrix[9, :] = bin9[1, :]
    #binMatrix[10, :] = bin10[1, :]

    #binMatrix = BestCutsOnTree(trees[2], hypergraph, incidence_struct, fixed_vertices, capacities, bestsolns_on_tree)
    #binMatrix = BestCutsOnTree(tree2, hypergraph, incidence_struct, fixed_vertices, capacities, bestsolns_on_tree)
    
    return binMatrix=#

    #=
    (m, ~) = size(binMatrix)
    union_cut = Int[]
    intersect_cut = Vector{Int}(1:hypergraph.e)

    for i in 1:m
        bins_i = binMatrix[i, :]
        (cut_edges_mrk, ~, ~, ~, ~) = cutProfile(hypergraph, incidence_struct, bins_i)
        cut_edges_i = findall(!iszero, cut_edges_mrk)
        union!(union_cut, cut_edges_i)
        intersect!(intersect_cut, cut_edges_i)
    end

    hsizes = hypergraph.eptr[2:end] - hypergraph.eptr[1:end-1]
    threshold_hedges = findall(x-> x > 3, hsizes)
    union!(union_cut, threshold_hedges)

    (cc, cs) = hypergraphCC(hypergraph, union_cut, false)
    vwts_cc = ContractVtxWts(hypergraph.vwts, cc)
    n_cc, e_cc, hedges_cc, eptr_cc, hwts_cc = ContractHyperGraph(hypergraph, cc)
    hypergraph_cc = Hypergraph(n_cc, e_cc, hedges_cc, eptr_cc, vwts_cc, hwts_cc)

    (adj_mat_cc, ~) = HypergraphToGraph(hypergraph_cc, expander_cycles)
    incidence_struct_cc = HypergraphToIncidence(hypergraph_cc)
    fixed_cc = -ones(Int, hypergraph_cc.n)

    for i in 1:length(fixed_part)
        if fixed_part[i] > -1
            fixed_cc[cc[i]] = fixed_part[i]
        end
    end

    p_0_cc = findall(x-> x == 0, fixed_cc)
    p_1_cc = findall(x-> x == 1, fixed_cc)
    fixed_vertices_cc = Pindex(p_0_cc, p_1_cc)
    X_cc = GenEigenVecs(hypergraph_cc, adj_mat_cc, fixed_vertices_cc, false, eigenvecs, eigen_iters)
    (tree_cc, ~) = GenEigenTree(adj_mat_cc, X_cc, eigenvecs, true)
    binMatrix_cc = BestCutsOnTree(tree_cc, hypergraph_cc, incidence_struct_cc, fixed_vertices_cc, capacities, 5)
    
    return binMatrix_cc
    =#
end

function SpectralCommunities(h_c::Hypergraph_C, incidence_struct::Incidence, ubfactor::Int; eigenvecs::Int = 2, expander_cycles::Int = 3, eigen_iters::Int = 20, bestsolns_on_tree::Int = 10)
    @info "DETECTING SPECTRAL COMMUNITIES"
    line_log = repeat("=", 80)
    @info "$line_log"

    hypergraph = h_c.HG
    fixed_part = h_c.fixed_part
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    (adj_mat, ~) = HypergraphToGraph(hypergraph, expander_cycles)
    X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)

    #return X
    #k_means_token =  Clustering.kmeans(X', 100)

    k_means_token =  ParallelKMeans.kmeans(Yinyang(), X', 100)

    return k_means_token #.assignments
end
