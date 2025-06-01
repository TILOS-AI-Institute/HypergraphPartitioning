module SpectralRefinement

using Shuffle
using JuMP
#using Cbc
using Gurobi
using LinearAlgebra
using DataStructures
using SparseArrays
using Random
using Statistics
using LightGraphs
using SimpleWeightedGraphs

include("PartitionStructures.jl")
#include("SpectralCommunityDetection/SpectralCommunities.jl")
include("SpectralCommunityDetection/SpectralPartitioning.jl")
include("GenerateHyperedgesHash.jl")
include("HypergraphToPairList.jl")
include("StructureDetection.jl")
include("HypergraphToIncidence.jl")
include("ReadFiles.jl")
include("PartitionILP.jl")
include("Coarsen.jl")
include("Uncoarsen.jl")
include("GenerateCoarseHypergraph.jl")
include("ExportHypergraph.jl")
include("EstimateClusteringMistakes.jl")
include("InitialPartitions.jl")
include("FM.jl")
include("WriteClusters.jl")
#include("ConvertXilinxFormatToHmetis.jl")

function kahypar(fname::String)
    cmd = `./KaHyPar -h $fname -k 2 -e 0.1 -o cut -m recursive -p ../kahypar/config/cut_rKaHyPar_sea20.ini`
    run(cmd, wait=true)
end

function hmetis(fname::String, Nparts::Int, UBfactor::Int, Nruns::Int, CType::Int, RType::Int, Vcycle::Int, Reconst::Int, dbglvl::Int, hmetis_exec::String)
    hmetis_log = "hmetis_log.txt"
    hmetis_cmd = hmetis_exec * " " * fname * " " * string(Nparts) * " " * string(UBfactor) * " " * string(Nruns) * " " * string(CType) * " " * string(RType) * " " * string(Vcycle) * " " * string(Reconst) * " " * string(dbglvl) * " > " * hmetis_log
    #hmetis_cmd = `$hmetis_exec $fname $Nparts $UBfactor $Nruns $CType $RType $Vcycle $Reconst $dbglvl`
    run(`sh -c $hmetis_cmd`, wait=true)
    cmd = "rm " * hmetis_log
    run(`sh -c $cmd`, wait=true)
end

function Golden_Evaluator(hypergraph_file::String, pfile::String, ub_factor::Int)
    (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = ReadHypergraphFile(hypergraph_file)
    hypergraph = Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)
    pvec = zeros(Int, num_vertices)

    f = open(pfile, "r")
    vi = 0

    for ln in eachline(f)
        vi += 1
        pvec[vi] = parse(Int, ln)
    end

    total_vwts = sum(vertex_weights)
    max_capacity = Int(round(total_vwts/2 + total_vwts * ub_factor * 0.01))
    min_capacity = Int(round(total_vwts/2 - total_vwts * ub_factor * 0.01))
    original_cutsize = Main.SpectralRefinement.GoldenEvaluator(hypergraph, pvec, max_capacity, min_capacity)
end


function GoldenEvaluator(hypergraph::Hypergraph, partition::Vector{Int}, max_capacity::Int, min_capacity::Int)
    cut_size = 0
    balance = false
    area_part = zeros(Int, 2)

    for i in 1:hypergraph.e
        start_idx = hypergraph.eptr[i]
        end_idx = hypergraph.eptr[i+1]-1
        base_part = partition[hypergraph.hedges[start_idx]]
        
        for j in start_idx+1:end_idx
            part = partition[hypergraph.hedges[j]]

            if base_part != part
                cut_size += hypergraph.hwts[i]
                break
            end
        end
    end

    for i in 1:length(partition)
        area_part[partition[i]+1] += hypergraph.vwts[i]
    end

    if area_part[1] >= min_capacity && area_part[1] <= max_capacity && area_part[2] >= min_capacity && area_part[2] <= max_capacity
        balance = true
    end

    @info "GOLDEN DETAILS: $cut_size :: $balance"

    return cut_size
end

function OverlayBasedClusteringAndILP(hg::String, ub_factor::Int, seeds::Vector{Int})
    hg_split = split(hg, "/")
    hg_name = hg_split[end]
    hg_name_clustered = "clustered_" * hg_name
    (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = ReadHypergraphFile(hg)
    fixed_vtxs = -ones(Int, num_vertices)
    hypergraph = Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)
    (hypergraph_processed, original_indices, new_indices, unused_indices) = IsolateIslands(hypergraph)
    fixed_vtxs_processed = fixed_vtxs[original_indices] 
    incidence_struct = HypergraphToIncidence(hypergraph_processed)
    incidence_list = GenerateIncidenceList(incidence_struct)
    hyperedge_pair_list = GenerateHypergraphPairList(hypergraph_processed)
    hyperedges_hash = GenerateHyperedgesHash(hypergraph_processed)
    fixed_vertex_flag  = maximum(fixed_vtxs_processed) > -1 ? true : false
    community = ones(Int, hypergraph_processed.n)
    hypergraph_c = Hypergraph_C(hypergraph_processed, incidence_list, hyperedge_pair_list, community, fixed_vtxs_processed, fixed_vertex_flag, hyperedges_hash)
    partition_matrix = zeros(Int, 5, hypergraph_processed.n)
    cut_list = zeros(Int, 5)

    part_name = "default_hmetis_rpt_julia_titan_500_runs/" * hg_name * ".k.2.UBfactor.10.seed."

    for i in 1:length(seeds)
        j = 0
        partition_file = part_name * string(seeds[i])

        f = open(partition_file, "r")

        for ln in eachline(f)
            j += 1
            cc = new_indices[j]

            if cc == 0
                continue
            end

            partition_matrix[i, cc] = parse(Int, ln)
        end

        cut_list[i] = findCutsize(partition_matrix[i, :], hypergraph_processed, incidence_struct)
    end

    #@info "CUTS :: $(cut_list)"
    t_clus = @elapsed  (hypergraph_clustered, incidence_struct_clustered, fixed_part_clustered, ~, ~, ~) = GenerateSpectralCommunities(partition_matrix, hypergraph_c, incidence_struct)
    partition_vector = zeros(Int, hypergraph_clustered.n)
    ExportHypergraph(hypergraph_clustered, hg_name_clustered, fixed_part_clustered)
    #@info "RUNNING CPLEX AS GOLDEN PARTITIONER"

    k = 0
        
    if hypergraph_clustered.e <= 2000

        cmd = "zhiang_for_bodhi/ilp_cplex_k_way/build/ilp_k_solver" * " " * hg_name_clustered * " " * string(2) * " " * string(ub_factor) * " > ilp_log.txt"
        t_part = @elapsed begin
            run(`sh -c $cmd`, wait=true)
            pname = hg_name_clustered * ".part." * string(2)
            f = open(pname, "r")
            
            for ln in eachline(f)
                k += 1
                partition_vector[k]  = parse(Int, ln)
            end

            close(f)

            cmd = "rm " * pname
            run(`sh -c $cmd`, wait=true)
            cmd = "rm ilp_log.txt"
            run(`sh -c $cmd`, wait=true)
        end
    else
        #@info "RUNNING HMETIS AS GOLDEN PARTITIONER"
        t_part = @elapsed hmetis(hg_name_clustered, 2, ub_factor, 10, 1, 1, 0, 1, 0, "./hmetis") #"SpectralCommunityDetection/hmetis")
        pname = hg_name_clustered * ".part." * string(2)
        f = open(pname, "r")

        for ln in eachline(f)
            k += 1
            partition_vector[k]  = parse(Int, ln)
        end

        close(f)

        cmd = "rm " * pname
        run(`sh -c $cmd`, wait=true)
    end
    
    cmd = "rm " * hg_name_clustered
    run(`sh -c $cmd`, wait=true)
    cut_size = findCutsize(partition_vector, hypergraph_clustered, incidence_struct_clustered)

    return cut_size, t_clus+t_part
end

function SpectralHmetisRefinement(;refine_iters::Int = 4, solver_iters::Int = 20, hg::String = "", pfile::String = "", fg::String = "", Nparts::Int = 2, hyperedges_threshold::Int = 900, ub::Int = 5, nev::Int = 1, cycles::Int = 1, seed::Int = 0, best_solns::Int=10, pseed::Int=0)
    BLAS.set_num_threads(Threads.nthreads())
    hg_split = split(hg, "/")
    hg_name = hg_split[end]
    hypergraph_file = hg
    fixed_file = fg
    num_parts = Nparts
    ub_factor = ub
    Random.seed!(seed)
    
    line_log = repeat("=", 80)
    @info "$line_log"
    @info "STARTING SUPERVISED SPECTRAL PARTITIONING ENGINE"
    @info "$line_log"

    t_elapsed_io = @elapsed begin
        (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = ReadHypergraphFile(hypergraph_file)
        fixed_vtxs = -ones(Int, num_vertices)
        
        if isempty(fg) == false
            ReadHypergraphFixedFile(fixed_file, fixed_vtxs)
        end
        
        hypergraph = Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)
        hsizes = hypergraph.eptr[2:end] - hypergraph.eptr[1:end-1]
        hsize_max = maximum(hsizes)
        (hypergraph_processed, original_indices, new_indices, unused_indices) = IsolateIslands(hypergraph)
        fixed_vtxs_processed = fixed_vtxs[original_indices] 
        incidence_struct = HypergraphToIncidence(hypergraph_processed)
        incidence_list = GenerateIncidenceList(incidence_struct)
        hyperedge_pair_list = GenerateHypergraphPairList(hypergraph_processed)
        hyperedges_hash = GenerateHyperedgesHash(hypergraph_processed)
    end

    max_capacity = Int(ceil(sum(vertex_weights) * (50+ub_factor)/100))
    min_capacity = sum(vertex_weights) - max_capacity
    total_vwts = sum(vertex_weights)

    @info "TOTAL VERTICES: $num_vertices"
    @info "TOTAL HYPEREDGES: $num_hyperedges"
    @info "POST PROCESSING :: TOTAL VERTICES: $(hypergraph_processed.n)"
    @info "POST PROCESSING :: TOTAL HYPEREDGES: $(hypergraph_processed.e)"
    @info "SIZE OF LARGEST HYPEREDGE: $hsize_max"
    @info "MAX CAPACITY CONSTRAINT: $max_capacity"
    @info "MIN CAPACITY CONSTRAINT: $min_capacity"
    @info "$line_log"

    config = "cut_rKaHyPar_sea20.ini"
    fixed_vertex_flag  = maximum(fixed_vtxs_processed) > -1 ? true : false
    community = ones(Int, hypergraph_processed.n)
    hypergraph_c = Hypergraph_C(hypergraph_processed, incidence_list, hyperedge_pair_list, community, fixed_vtxs_processed, fixed_vertex_flag, hyperedges_hash)
    partition_vector = zeros(Int, hypergraph_processed.n)
    final_part = zeros(Int, num_vertices)

    t_hint = @elapsed begin 
        if isempty(pfile) == true
            i = 0

            @info "RUNNING HMETIS ON ORIGINAL HG TO GENERATE HINT FOR SPECTRAL"

            run(`hmetis $hypergraph_file "" $num_parts $ub_factor 10 1 1 0 1 0 $pseed`, wait=true)
            #run(`hmetis_standalone/hmetis $hypergraph_file $num_parts $ub_factor 10 1 1 0 1 0`, wait=true)

            #hmetis(hg, num_parts, ub_factor, 10, 1, 1, 0, 1, 0, "zhiang_for_bodhi/hmetis/hmetis")
            pname = hg * ".part." * string(num_parts)

            f = open(pname, "r")
            v = 0
            for ln in eachline(f)
                v += 1
                final_part[v] = parse(Int, ln)
                cc = new_indices[v]

                if cc == 0
                    continue
                end

                partition_vector[cc] = final_part[v]
            end

            close(f)

            #=f = open(pname, "r")
            for ln in eachline(f)
                part_i = split(ln)

                for j in 2:length(part_i)
                    v = parse(Int, part_i[j])
                    final_part[v] = i
                    cc = new_indices[v]

                    if cc == 0
                        continue
                    end

                    partition_vector[cc] = i
                end
                i += 1
            end
            
            close(f)=#

            cmd = "rm " * pname
            run(`sh -c $cmd`, wait=true)
        else
            i = 0
            f = open(pfile, "r")

            for ln in eachline(f)
                i += 1
                final_part[i] = parse(Int, ln)
                cc = new_indices[i]
        
                if cc == 0 
                    continue
                end
        
                partition_vector[cc]  = parse(Int, ln)
            end

            close(f)
        end
    end

    hmetis_cut = findCutsize(partition_vector, hypergraph_processed, incidence_struct)

    partition_vector_orig = deepcopy(partition_vector)

    @info "[HINT] CUT RECORDED IS $hmetis_cut"

    f = open("twgts.wts", "w")
    t_ub = (50 - ub_factor)/100
    println(f, "0 = " * string(t_ub))
    println(f, "1 = " * string(1.0-t_ub))
    close(f)

    t_spec = @elapsed best_part_matrix, best_cuts = RefineIteratively(partition_vector, hypergraph_c, incidence_struct, ub_factor, expander_cycles = cycles, eigenvecs = nev, iters=refine_iters, eigen_iters = solver_iters, mseed=seed)

    @info "Cuts in pool: $best_cuts"

    m, n = size(best_part_matrix)
    best_solns = best_solns > m ? m : best_solns
    cuts_perm = sortperm(best_cuts)
    best_cuts = best_cuts[cuts_perm[1:best_solns]]
    best_part_matrix = best_part_matrix[cuts_perm[1:best_solns], :]

    @info "Cuts picked: $best_cuts"

    global_cut, global_tree = findmin(best_cuts)
    global_part = best_part_matrix[global_tree, :]
    

    t_clus = @elapsed (hypergraph_clustered, incidence_struct_clustered, fixed_part_clustered, vwts_clustered, cclist, csize) = GenerateSpectralCommunities(best_part_matrix, hypergraph_c, incidence_struct)

    @info "$line_log"
    @info "SIZE OF CLUSTERED HYPERGRAPH IS: $(hypergraph_clustered.n) VERTICES AND $(hypergraph_clustered.e) HYPEREDGES"
    @info "$line_log"
    @info "PARTITIONING CLUSTERED HYPERGRAPH ....."

    fixed_part_clustered = -ones(Int, hypergraph_clustered.n)
    capacities = [max_capacity, max_capacity]
    partition_vector = zeros(Int, hypergraph_clustered.n)
    t_part = 0
    hg_name_clustered = "clustered_" * hg_name
    ExportHypergraph(hypergraph_clustered, hg_name_clustered, fixed_part_clustered)
    i = 0

    if hypergraph_clustered.e < hyperedges_threshold
        @info "RUNNING CPLEX AS GOLDEN PARTITIONER"
        cmd = "zhiang_for_bodhi/ilp_cplex_k_way/build/ilp_k_solver" * " " * hg_name_clustered * " " * string(num_parts) * " " * string(ub_factor) * " > ilp_log.txt"
        t_part = @elapsed begin
            run(`sh -c $cmd`, wait=true)
            pname = hg_name_clustered * ".part." * string(Nparts)
            f = open(pname, "r")
            
            for ln in eachline(f)
                i += 1
                partition_vector[i]  = parse(Int, ln)
            end

            close(f)

            cmd = "rm " * pname
            run(`sh -c $cmd`, wait=true)
            cmd = "rm ilp_log.txt"
            run(`sh -c $cmd`, wait=true)
        end
    else
        @info "RUNNING HMETIS AS GOLDEN PARTITIONER"
        t_part = @elapsed hmetis(hg_name_clustered, num_parts, ub_factor, 10, 1, 1, 0, 1, 0, "./SpecPart/hmetis") #"SpectralCommunityDetection/hmetis")
        pname = hg_name_clustered * ".part." * string(Nparts)
        f = open(pname, "r")

        for ln in eachline(f)
            i += 1
            partition_vector[i]  = parse(Int, ln)
        end

        close(f)

        cmd = "rm " * pname
        run(`sh -c $cmd`, wait=true)
    end

    cmd = "rm " * hg_name_clustered
    run(`sh -c $cmd`, wait=true)
    
    part_area = [0, 0]
    partition_vector_top = zeros(Int, hypergraph_processed.n)
    post_tool_cut = findCutsize(partition_vector, hypergraph_clustered, incidence_struct_clustered)

    @info "[POST TOOL CUT] CUT RECORDED IS $post_tool_cut"

    tool_cut = false

    if post_tool_cut < global_cut
        tool_cut = true
        global_cut = post_tool_cut
        global_part = partition_vector
    end
        
    if global_cut >= hmetis_cut

        f = open(hg_name * "_" * string(ub_factor) * ".part." * string(num_parts), "w")

        for i in 1:length(final_part)
            p = final_part[i]

            part_area[p+1] += hypergraph.vwts[i]
            println(f, final_part[i])
        end

        close(f)

        global_cut = hmetis_cut
    else
        #GoldenEvaluator(hypergraph_processed, partition_vector_orig, max_capacity, min_capacity)
        #final_part = partition

        if tool_cut == true
            for i in 1:length(cclist)
                #part_area[global_part[cclist[i]]+1] += hypergraph_processed.vwts[i]
                partition_vector_top[i] = global_part[cclist[i]]
            end
        else
            partition_vector_top = global_part
        end 
        
        for i in 1:length(original_indices)  
            v = original_indices[i]
            final_part[v] = partition_vector_top[i]
            part_area[final_part[v]+1] += hypergraph.vwts[v]
        end

        for i in 1:length(unused_indices)
            v = unused_indices[i]

            min_side = part_area[1] < part_area[2] ? 1 : 2
            final_part[v] = min_side-1
            part_area[min_side] += vertex_weights[v]
        end

        f = open(hg_name * "_" * string(ub_factor) * ".part." * string(num_parts), "w")

        for i in 1:length(final_part)
            println(f, final_part[i])
        end

        close(f)

        #x = GoldenEvaluator(hypergraph, final_part, max_capacity, min_capacity)
    end

    @info "$line_log"
    @info "[SUPERVISED SPECTRAL] CUTSIZE OBTAINED: $global_cut"
    @info "[SUPERVISED SPECTRAL] AREA SPLIT OBTAINED: $(part_area[1]) [$(round(part_area[1]/total_vwts, digits=2))%] :: $(part_area[2]) [$(round(part_area[2]/total_vwts, digits=2))%]"
    @info "$line_log"
    @info "[RUNTIME] IO PROCESSING :: $t_elapsed_io SECONDS"
    @info "[RUNTIME] SPECTRAL :: $t_spec SECONDS"
    @info "[RUNTIME] CLUSTERING :: $t_clus SECONDS"
    @info "[RUNTIME] PARTITION CLUSTERED HYPERGRAPH :: $t_part SECONDS"
    @info "$line_log"
    @info "[RUNTIME] TOTAL EXECUTION TIME :: $(t_elapsed_io + t_spec + t_clus + t_part) SECONDS"
    @info "$line_log"

    return global_cut, hmetis_cut
end

end #Module ends here
