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

design_list = ["ibm01",
               "ibm02",
               "ibm03",
               "ibm04",
               "ibm05",
               "ibm06",
               "ibm07",
               "ibm08",
               "ibm09",
               "ibm10",
               "ibm11",
               "ibm12",
               "ibm13",
               "ibm14",
               "ibm15",
               "ibm16",
               "ibm17",
               "ibm18"
                ]

design_list = ["ibm01.weight",
                "ibm02.weight",
                "ibm03.weight",
                "ibm04.weight",
                "ibm05.weight",
                "ibm06.weight",
                "ibm07.weight",
                "ibm08.weight",
                "ibm09.weight",
                "ibm10.weight",
                "ibm11.weight",
                "ibm12.weight",
                "ibm13.weight",
                "ibm14.weight",
                "ibm15.weight",
                "ibm16.weight",
                "ibm17.weight",
                "ibm18.weight"
                 ]

num_parts_list = [2, 3, 4]
ub_factor_list = [2]
benchmark_dir = "/home/fetzfs_projects/SpecPart/testcases/ibm_w"

for design in design_list
    for num_parts in num_parts_list
        regression_file_name = "/home/bodhi91/SpecPart/regression/" * string(num_parts) * "_way/regression_results.wts.csv" 
        f = open(regression_file_name, "w")
            println(f, "Testcase, Numparts, UBfactor, Cutsize, Runtime")
        close(f)
    end
end


for num_parts in num_parts_list
    regression_file_name = "/home/bodhi91/SpecPart/regression/" * string(num_parts) * "_way/regression_results.wts.csv" 
    hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
    for ub_factor in ub_factor_list
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
