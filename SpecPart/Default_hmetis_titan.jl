include("SpectralRefinement.jl")

partitioner = "hmetis"
benchmark = "Titan"

UBfactor_list = [10]
Nparts_list = [2]
eig_vecs_list = [2] # number of eigenvalues
num_cycles_list = [2]  # number of cycles used for converting hypergraph
h_threshold_list = [300] # threshold of number of hyperedges in the clustered_hypergraph for ILP solver
solver_iters_list = [80] # number of iterations in the eigenvalue solvers
#spec_iters_list = [3, 2, 1] # number of iterations for the spectral framework
#best_solns_list = [10, 8, 6, 4, 2] # number of best partition solutions used for cut elimination
spec_iters_list = [2] #, 3, 4, 5] # number of iterations for the spectral framework
best_solns_list = [5] # number of best partition solutions used for cut elimination
seed_list = Vector{Int}(1:499)

work_dir = pwd()
work_dir *= "/default_hmetis_rpt_julia_titan_500_runs/"

if isdir(work_dir) == false
    mkdir(work_dir)
end

result_file = work_dir * partitioner * "_" * benchmark * "_summary_file.txt"
julia_exe = "/home/bodhi91/sandbox/software/julia-1.6.0-linux-x86_64/julia-1.6.0/bin/julia"
hmetis_exe = pwd() * "/zhiang_for_bodhi/hmetis/hmetis"
patoh_exe = pwd() * "/patoh"
kahypar_exe = pwd() * "/KaHyPar"
kahypar_config = "/home/bodhi91/sandbox/Partitioners/IterativeSpectralRefinement/cut_rKaHyPar_sea20.ini"
design_list = String[]
benchmark_dir = "/home/bodhi91/sandbox/Partitioners/Benchmarks/" * benchmark * "/"

if benchmark == "IBM"
    design_list = ["" for i in 1:18]
    for i in 1:length(design_list)
        if i < 10
            design_list[i] = "ibm0" * string(i)
        else
            design_list[i] = "ibm"  * string(i)
        end
    end
elseif benchmark == "IBM_w"
    design_list = ["" for i in 1:18]
    for i in 1:length(design_list)
        if i < 10
            design_list[i] = "ibm0" * string(i) * ".weight"
        else
            design_list[i] = "ibm"  * string(i) * ".weight"
        end
    end
elseif benchmark == "Titan"
    design_list = ["gsm_switch", "sparcT2_core"]
    #=
    design_list = ["bitcoin_miner",
                "bitonic_mesh",
                "cholesky_bdti",
                "cholesky_mc",
                "dart",
                "denoise",
                "des90",
                "directrf",
                "gsm_switch",
                "LU230",
                "LU_Network",
                "mes_noc",
                "minres",
                "neuron",
                "openCV",
                "segmentation",
                "SLAM_spheric",
                "sparcT1_chip2",
                "sparcT1_core",
                "sparcT2_core",
                "stap_qrd",
                "stereo_vision"]=#
else
    print("[INFO] [ERROR] Invalid Benchmarks!!!")
end

f = open(result_file, "w")
line  = "design,Nparts,UBfactor,Seed,"
line *= "cutsize,runtime\n"
println(f, line)
close(f)

function ConvertPartitionFileFormat(pfile::String)
    f = open(pfile, "r")
    n = 0
    ivec_arr = []

    for ln in eachline(f)
        ivec = split(ln)
        push!(ivec_arr, ivec[2:end])
        n += length(ivec) - 1
    end
    close(f)

    pvec = zeros(Int, n)

    for i in length(ivec_arr)
        vtxs = parse.(Int, ivec_arr[i])
        pvec[vtxs] .= i-1
    end

    cmd = "rm " * pfile
    run(`sh -c $cmd`, wait=true)

    f = open(pfile, "w")

    for p in pvec
        println(f, p)
    end

    close(f)
end

for design in design_list
    for Nparts in Nparts_list
        for UBfactor in UBfactor_list
            for seed in seed_list
                original_hypergraph_file = benchmark_dir * design * ".hgr"
                hypergraph_file = work_dir * design * ".hgr"
                cmd = "cp " * original_hypergraph_file * " " * hypergraph_file
                run(`sh -c $cmd`, wait=true)
                t = @elapsed run(`$hmetis_exe $hypergraph_file "" $Nparts $UBfactor 10 1 1 0 1 0 $seed`, wait=true)
                #ConvertPartitionFileFormat(hypergraph_file * ".part." * string(Nparts))
                original_runtime = round(t, digits=2)
                original_initial_solution_file = hypergraph_file * ".part." * string(Nparts)
                initial_solution_file = hypergraph_file * ".k." * string(Nparts) *  ".UBfactor." * string(UBfactor) * ".seed." * string(seed)
                cmd = "mv  " * original_initial_solution_file * " " * initial_solution_file
                run(`sh -c $cmd`, wait=true)
                (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = Main.SpectralRefinement.ReadHypergraphFile(hypergraph_file)
                pvec = zeros(Int, num_vertices)
                p = open(initial_solution_file, "r")
                vi = 0

                for ln in eachline(p)
                    vi += 1
                    pvec[vi] = parse(Int, ln)
                end

                close(p)

                total_vwts = sum(vertex_weights)
                max_capacity = Int(round(total_vwts/Nparts + total_vwts * UBfactor * 0.01))
                min_capacity = Int(round(total_vwts/Nparts - total_vwts * UBfactor * 0.01))
                hypergraph = Main.SpectralRefinement.Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)
                original_cutsize = Main.SpectralRefinement.GoldenEvaluator(hypergraph, pvec, max_capacity, min_capacity)
    
                println("[INFO]  Design : ",  design, " UBfactor : ", UBfactor, "   Nparts : ", Nparts)
                println("[INFO]  Original Summary   CutSize : ", original_cutsize, "    Runtime : ", original_runtime)

                line = design * "," * string(Nparts)
                line *= "," * string(UBfactor)
                line *= "," * string(seed)
                line *= "," * string(original_cutsize)
                line *= "," * string(original_runtime)
                line *= "\n"
                f = open(result_file, "a")
                write(f, line)
                close(f)
            end
        end
    end
end
