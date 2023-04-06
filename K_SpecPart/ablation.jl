include("specpart.jl")

benchmark_dir = "/home/fetzfs_projects/SpecPart/testcases/titan"

design_list = ["sparcT1_core",
               "cholesky_mc",
               "segmentation",
               "denoise",
               "gsm_switch"]

hmetis_cuts = [983 1890 2569 2892 ; 281 546 979 1330 ; 85 376 487 598 ; 441 820 1040 1137 ; 3786 5255 3930 5500 ;]

num_parts_list = [2, 3, 4, 5]
ub_factor = 5



# default values
eigenvecs_def = 2
best_solns_def = 5
solve_iters_def = 40
refine_iters_def = 2
cycles_def = 1

# effect of eigenvectors 

eigenvectors_list = [2, 3, 4]
best_solns_list = [5, 6, 7, 8]
solver_iters_list = [20, 40, 60, 80]
refine_iters_list = [1, 2, 3, 4]
cycles_list = [1, 2, 3, 4]


# effect of eigenvectors 
function run_ablation_eigenvecs() 
    f = open("ablation_eigenvecs.dat", "w")
    for eigen_vecs in eigenvectors_list
        improvement = 0.0
        iters = 0
        for i in 1:length(design_list)
            design = design_list[i]
            hypergraph_file = benchmark_dir * "/" * design * ".hgr"                               
            for j in 1:length(num_parts_list) 
                num_parts = num_parts_list[j]
                hmetis_cut = hmetis_cuts[i, j]
                iters += 1
                hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
                hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
                (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                    hint_file=hint_file, 
                                                    imb = ub_factor, 
                                                    num_parts = num_parts, 
                                                    eigvecs = eigen_vecs, 
                                                    best_solns = best_solns_def, 
                                                    solver_iters = solve_iters_def, 
                                                    refine_iters = refine_iters_def,
                                                    ncycles = cycles_def)
                diff = (hmetis_cut - cut_token[1])/hmetis_cut * 100
                if diff < 0 
                    diff = 100.0 + (diff*-1.0)
                end
                improvement += diff
            end
        end
        improvement /= iters
        line = string(eigen_vecs) * "," * string(improvement)
        println(f, line) 
    end
    close(f)
end

# effect of best solns
function run_ablation_best_solns() 
    f = open("ablation_best_solns.dat", "w")
    for best_solns in best_solns_list
        improvement = 0.0
        iter = 0
        for i in 1:length(design_list)
            design = design_list[i]
            hypergraph_file = benchmark_dir * "/" * design * ".hgr"                               
            for j in 1:length(num_parts_list) 
                iters += 1
                num_parts = num_parts_list[j]
                hmetis_cut = hmetis_cuts[i, j]
                hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
                hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
                (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                    hint_file=hint_file, 
                                                    imb = ub_factor, 
                                                    num_parts = num_parts, 
                                                    eigvecs = eigen_vecs_def, 
                                                    best_solns = best_solns, 
                                                    solver_iters = solve_iters_def, 
                                                    refine_iters = refine_iters_def,
                                                    ncycles = cycles_def)
                diff = (hmetis_cut - cut_token[1])/hmetis_cut * 100
                if diff < 0 
                    diff = 100.0 + (diff*-1.0)
                end
                improvement += diff
            end
        end
        improvement /= iters
        line = string(best_solns) * "," * string(improvement)
        println(f, line) 
    end
    close(f)
end

# effect of solver iters
function run_ablation_solver_iters() 
    f = open("ablation_solver_iters.dat", "w")
    for solver_iters in solver_iters_list
        improvement = 0.0
        iters = 0
        for i in 1:length(design_list)
            design = design_list[i]
            hypergraph_file = benchmark_dir * "/" * design * ".hgr"                               
            for j in 1:length(num_parts_list) 
                iters += 1
                num_parts = num_parts_list[j]
                hmetis_cut = hmetis_cuts[i, j]
                hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
                hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
                (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                    hint_file=hint_file, 
                                                    imb = ub_factor, 
                                                    num_parts = num_parts, 
                                                    eigvecs = eigen_vecs_def, 
                                                    best_solns = best_solns_def, 
                                                    solver_iters = solver_iters, 
                                                    refine_iters = refine_iters_def,
                                                    ncycles = cycles_def)
                diff = (hmetis_cut - cut_token[1])/hmetis_cut * 100
                if diff < 0 
                    diff = 100.0 + (diff*-1.0)
                end
                improvement += diff
            end
        end
        improvement /= iters
        line = string(eigen_vecs) * "," * string(improvement)
        println(f, line)
    end
    close(f)
end

# effect of specpart iters
function run_ablation_specpart_iters() 
    f = open("ablation_specpart_iters.dat", "w")
    for refine_iters in refine_iters_list
        improvement = 0.0
        iters = 0
        for i in 1:length(design_list)
            design = design_list[i]
            hypergraph_file = benchmark_dir * "/" * design * ".hgr"                               
            for j in 1:length(num_parts_list) 
                iters += 1
                num_parts = num_parts_list[j]
                hmetis_cut = hmetis_cuts[i, j]
                hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
                hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
                (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                    hint_file=hint_file, 
                                                    imb = ub_factor, 
                                                    num_parts = num_parts, 
                                                    eigvecs = eigen_vecs_def, 
                                                    best_solns = best_solns_def, 
                                                    solver_iters = solver_iters_def, 
                                                    refine_iters = refine_iters,
                                                    ncycles = cycles_def)
                diff = (hmetis_cut - cut_token[1])/hmetis_cut * 100
                if diff < 0 
                    diff = 100.0 + (diff*-1.0)
                end
                improvement += diff
            end
        end
        improvement /= iters
        line = string(solver_iters) * "," * string(improvement)
        println(f, line)
    end
    close(f)
end

# effect of cycles 
function run_ablation_cycles() 
    f = open("ablation_cycles.dat", "w")
    for cycles in cycles_list
        improvement = 0.0
        iters = 0
        for i in 1:length(design_list)
            design = design_list[i]
            hypergraph_file = benchmark_dir * "/" * design * ".hgr"                               
            for j in 1:length(num_parts_list) 
                iters += 1
                num_parts = num_parts_list[j]
                hmetis_cut = hmetis_cuts[i, j]
                hint_dir = "/home/fetzfs_projects/SpecPart/specpart_hints/" * string(num_parts) * "_way"
                hint_file = hint_dir * "/ub_factor_" * string(ub_factor) * "/" * design * ".hgr.part." * string(num_parts)
                (~, cut_token) = Main.SpecPart.specpart_run(hypergraph_file, 
                                                    hint_file=hint_file, 
                                                    imb = ub_factor, 
                                                    num_parts = num_parts, 
                                                    eigvecs = eigen_vecs_def, 
                                                    best_solns = best_solns_def, 
                                                    solver_iters = solver_iters_def, 
                                                    refine_iters = refine_iters_def,
                                                    ncycles = cycles)
                diff = (hmetis_cut - cut_token[1])/hmetis_cut * 100
                if diff < 0 
                    diff = 100.0 + (diff*-1.0)
                end
                improvement += diff
            end
        end
        improvement /= iters
        line = string(cycles) * "," * string(improvement)
        println(f, line)
    end
    close(f)
end

run_ablation_eigenvecs()
run_ablation_best_solns()
run_ablation_cycles()
run_ablation_solver_iters()
run_ablation_specpart_iters()