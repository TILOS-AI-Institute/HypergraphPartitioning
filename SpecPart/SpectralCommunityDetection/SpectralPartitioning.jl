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
using Match
using Dates
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

function metis(fname::String, seed::Int, opts::String)
    time_stamp = string(now())
    log_name = "metis" * opts * "." * time_stamp * "log.txt"
    twname = "twgts.wts"
    #metis_r = "./gpmetis " * fname * " 2 -ptype=rb -ufactor=100 -dbglvl=0" # > 'metis_log.txt'"""
    #metis_r = `./gpmetis $fname 2 -ptype=rb -ufactor=100 -dbglvl=0` 
    metis_r = `./metis_script.sh $fname $log_name $twname $seed`
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

function ComputeTreePartition(adj_mat::SparseMatrixCSC, X::Array{Float64}, wt_matrix::SparseMatrixCSC, hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vertices::Pindex, capacities::Vector{Int}, opts::Int, seed::Int)
    @match opts begin
        1 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 1], 2)
        end

        2 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 1], 1, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 1], 1)
        end

        3 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 1], 1, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 1], 4)
        end

        4 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 2], 2)
        end

        5 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 2], 1, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 2], 1)
        end

        6 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 2], 1, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 2], 4)
        end

        7 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X, 2, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X, 2)
        end

        8 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X, 2, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X, 1)
        end

        9 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 3], 1, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 3], 2)
        end

        10 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 3], 1, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 3], 1)
        end

        11 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, 3], 1, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, 3], 4)
        end

        12 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, [1,3]], 2, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, [1,3]], 2)
        end

        13 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, [1,3]], 2, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, [1,3]], 1)
        end

        14 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, [2,3]], 2, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, [2,3]], 2)
        end

        15 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X[:, [2,3]], 2, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X[:, [2,3]], 1)
        end

        16 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X, 3, false)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X, 2)
        end

        17 => begin
        clique_graph = ModifyExpanderWts(adj_mat, X, 3, true)
        tree, treeMatrix = GenEigenTree(clique_graph, wt_matrix, X, 1)
        end

        _ => @warn "Enter Correct Tree Opts!!"
    end

    pvec_metis = MetisOnTree(hypergraph, incidence_struct, fixed_vertices, tree, treeMatrix, capacities, seed, opts)
    cut_metis = findCutsize(pvec_metis, hypergraph, incidence_struct)
    pvec_sweep, cut_sweep = FindBestCutOnTree(tree, hypergraph, incidence_struct, fixed_vertices, capacities, 1)

    #@info "Tree type: $opts :: cut recorded :: $cut_metis :: $cut_sweep"

    if cut_metis < cut_sweep
        return pvec_metis, cut_metis
    else
        return pvec_sweep, cut_sweep
    end
end

function METIS_Run(hypergraph::Hypergraph, tree::SimpleWeightedGraph, seed::Int, opts::String)
    fname = BuildMetisGraph(tree, opts)
    metis(fname, seed, opts)

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

function MetisOnTree(hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vertices::Pindex, tree::SimpleWeightedGraphs.SimpleGraph, treeMatrix::SparseMatrixCSC, capacities::Vector{Int}, seed::Int, opts::Int)
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

    cut_elements, pvec = METIS_Run(hypergraph, g, seed, string(opts))

    #@info "METIS CUT FOUND IS: $(cut_elements[1]) WITH SPLIT: $(cut_elements[2])"

    return pvec
end

function RefineIteratively(pvec::Vector{Int}, hypergraph_c::Hypergraph_C, incidence_struct::Incidence, ubfactor::Int; eigenvecs::Int = 1, expander_cycles::Int = 1, eigen_iters::Int = 50, iters::Int=10, mseed::Int=0)
    families = 0
    
    if eigenvecs == 1 
        families = 3
    elseif eigenvecs == 2
        families = 8
    else
        families = 17
    end

    hypergraph = hypergraph_c.HG
    min_cut = 1e10
    min_tree = 0
    all_tree_parts = zeros(Int, families*iters, hypergraph.n)
    tree_part_family = zeros(Int, families, hypergraph.n)
    all_cuts = zeros(Int, families*iters)
    list_of_cuts = zeros(Int, families)
    best_family = tree_part_family
    best_list_of_cuts = list_of_cuts
    global_min_cut = min_cut
    max_capacity = Int(ceil(sum(hypergraph.vwts) * (50+ubfactor)/100))
    capacities = [max_capacity, max_capacity]
    (adj_mat, wt_matrix) = HypergraphToGraph(hypergraph, 1)
    (adj_mat_tree, wt_matrix_tree) = HypergraphToGraph(hypergraph, expander_cycles)
    #=p_0 = findall(x-> x == 0, pvec)
    p_1 = findall(x-> x == 1, pvec)
    fixed_vertices = Pindex(p_0, p_1)
    X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)

    return X=#
    k = 0

    for i in 1:iters
        p_0 = findall(x-> x == 0, pvec)
        p_1 = findall(x-> x == 1, pvec)
        fixed_vertices = Pindex(p_0, p_1)
        X = GenEigenVecs(hypergraph, adj_mat, fixed_vertices, false, eigenvecs, eigen_iters)
        fixed_vertices = Pindex(Int[], Int[])
        #added

        t_iter = @elapsed Threads.@threads for j in 1:families
            tree_part_family[j, :], list_of_cuts[j] = ComputeTreePartition(adj_mat_tree, X, wt_matrix_tree, hypergraph, incidence_struct, fixed_vertices, capacities, j, mseed)  
        end

        min_cut, min_tree = findmin(list_of_cuts)
        @info "[SPECTRAL ITERATION $i] MIN CUT FOUND :: $min_cut :: TREE GEN-SOLVE TIME :: $t_iter seconds"
        
        for j in 1:families
            k += 1
            all_tree_parts[k, :] = tree_part_family[j, :]
            all_cuts[k] = list_of_cuts[j]
        end

        #=if min_cut <= global_min_cut
            best_family = deepcopy(tree_part_family)
            best_list_of_cuts = deepcopy(list_of_cuts)
            global_min_cut = min_cut
        en=#
        
        pvec = tree_part_family[min_tree, :]
    end

    #@info "[BEST SPECTRAL CUT :: $global_min_cut]"

    #cperm = sortperm(all_cuts)

    #return (all_tree_parts, all_cuts)

    return (best_family, best_list_of_cuts)
end
