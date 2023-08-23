include("specpart.jl")

design_list = ["sparcT1_core",
               "neuron",
               "stereo_vision",
               "des90",
               "SLAM_spheric",
               "cholesky_mc",
               "segmentation",
               "bitonic_mesh",
               "dart",
               "openCV",
               "stap_qrd",
               "minres",
               "cholesky_bdti",
               "denoise",
               "sparcT2_core",
               "gsm_switch",
               "mes_noc",
               "LU230",
               "LU_Network",
               "sparcT1_chip2",
               "directrf",
               "bitcoin_miner"]

num_parts_list = [3, 4]
ub_factor_list = [2]
benchmark_dir = "/home/fetzfs_projects/SpecPart/testcases/titan"

for design in design_list
    for num_parts in num_parts_list
        regression_file_name = "/home/fetzfs_projects/SpecPart/regression/" * string(num_parts) * "_way/regression_results.csv" 
        f = open(regression_file_name, "w")
            println(f, "Testcase, Numparts, UBfactor, Cutsize, Runtime")
        close(f)
    end
end


for num_parts in num_parts_list
    regression_file_name = "/home/fetzfs_projects/SpecPart/regression/" * string(num_parts) * "_way/regression_results.csv" 
    hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
    for ub_factor in ub_factor_list
        #=if num_parts == 2 && ub_factor == 20
            continue
        end
        if num_parts > 2 && ub_factor == 5
            continue
        end=#
        if num_parts == 3 && ub_factor != 20
            continue
        end
        for design in design_list
            hypergraph_file = benchmark_dir * "/" * design * ".hgr"
            hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
            t = @elapsed (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                hint_file=hint_file, 
                                                imb = ub_factor, 
                                                num_parts = num_parts, 
                                                eigvecs = 2, 
                                                best_solns = 5, 
                                                solver_iters = 40, 
                                                refine_iters = 2)
            cutsize = cut_token[1]
            @info "Cutsize recorded at end of specpart $cutsize" 
            f = open(regression_file_name, "a")
            line = design * ", " 
            line *= string(num_parts) * ", "
            line *= string(ub_factor) * ", "
            line *= string(cutsize) * ", " 
            line *= string(t)
            println(f, line)
            close(f)
        end
    end
end

#=

for design in design_list
    hypergraph_file = benchmark_dir * "/" * design * ".hgr"
    for num_parts in num_parts_list
        hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
        regression_file_name = "/home/fetzfs_projects/SpecPart/regression/" * string(num_parts) * "_way/regression_results.csv" 
        for ub_factor in ub_factor_list
            if num_parts == 2 && ub_factor == 20
                continue
            end
            if num_parts > 2 && ub_factor == 5
                continue
            end
            hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
            t = @elapsed (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                hint_file=hint_file, 
                                                imb = ub_factor, 
                                                num_parts = num_parts, 
                                                eigvecs = 2, 
                                                best_solns = 5, 
                                                solver_iters = 40, 
                                                refine_iters = 2)
            cutsize = cut_token[1]
            @info "Cutsize recorded at end of specpart $cutsize" 
            f = open(regression_file_name, "a")
            line = design * ", " 
            line *= string(num_parts) * ", "
            line *= string(ub_factor) * ", "
            line *= string(cutsize) * ", " 
            line *= string(t)
            println(f, line)
            close(f)
        end
    end
end
=#
