include("SpectralRefinement.jl")

hypergraph_files = ["../../SmartPartAndExRefine/titan23_benchmark/CUHK_benchmark/sparcT2_core.hgr"]
ub_factor = 2

starts = Vector{Int}(1:20);
seeds = Vector{Int}(0:19);
hmetis_exe = pwd() * "/zhiang_for_bodhi/hmetis/hmetis"

for i in 1:length(hypergraph_files)
    hypergraph_file = hypergraph_files[i]
    result_file = pwd() * "/hMetisVsSpecPart_" * string(i) * ".csv"
    f = open(result_file, "w")
    println(f, "#Starts,hMetis,SpecPart")
    close(f)
    
    for start in starts
        min_cutsize = 1e09

        for j in 1:start
            #Main.SpectralRefinement.hmetis(hypergraph_file, 2, ub_factor, 10, 1, 1, 0, 1, 0, "./hmetis")
            run(`$hmetis_exe $hypergraph_file "" 2 $ub_factor 10 1 1 0 1 0 $(seeds[j])`, wait=true)
            (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = Main.SpectralRefinement.ReadHypergraphFile(hypergraph_file)
            max_capacity = Int(ceil(sum(vertex_weights) * (50+ub_factor)/100))
            min_capacity = sum(vertex_weights) - max_capacity
            hypergraph = Main.SpectralRefinement.Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)  
            pvec = zeros(Int, num_vertices)
            original_part_file = hypergraph_file * ".part.2"
            fp = open(original_part_file, "r")
            k = 0

            for ln in eachline(fp)
                k += 1
                pvec[k] = parse(Int, ln)
            end

            close(fp)

            hmetis_cutsize = Main.SpectralRefinement.GoldenEvaluator(hypergraph, pvec, max_capacity, min_capacity)

            if hmetis_cutsize < min_cutsize
                min_cutsize = hmetis_cutsize
                cmd = "mv " * original_part_file * " " * "best_part.part.2"
                run(`sh -c $cmd`, wait=true)
            end
        end

        (specpart_cutsize, ~) = Main.SpectralRefinement.SpectralHmetisRefinement(hg = hypergraph_file, pfile = "best_part.part.2", cycles=2, hyperedges_threshold=600, ub=ub_factor, nev=2, seed=0, solver_iters=80, refine_iters=2, best_solns=3);
        
        f = open(result_file, "a")
        write(f, string(start), ",", string(min_cutsize), ",", string(specpart_cutsize), "\n")
        close(f)
    end
end


