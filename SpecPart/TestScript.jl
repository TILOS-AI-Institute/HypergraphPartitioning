using DataFrames, CSV

include("SpectralRefinement.jl")

#=files = readdir("/home/bodhi91/DEV/sandbox/Partitioner/Benchmarks_Spectral/titan23_benchmark/CUHK_benchmark/", join=true)
file_names = readdir("/home/bodhi91/DEV/sandbox/Partitioner/Benchmarks_Spectral/titan23_benchmark/CUHK_benchmark/")
nparts_list = [2]
ub_list = [10]
nev_list = [3]

#files = readdir("../SmartPartAndExRefine/ISPD_benchmark/hmetis_format/", join=true)

t = @elapsed for ub_fac in ub_list
    for n_eig in nev_list
        csv_name = "TITAN_" * string(ub_fac) * "_" * string(n_eig) * "_results.csv" 
        test_case_list = String[]
        spec_cut_list = Int[]
        hmetis_cut_list = Int[]

        for i in 1:length(files)
            partition_vector = Int[]
            f = files[i]
            push!(test_case_list, file_names[i])
            @info "Running $f with ub_fac = $ub_fac, eigen_vecs = $n_eig"
            ~, spec_cut, hmetis_cut = Main.SpectralRefinement.SpectralHmetisRefinement(hg = f, nparts = 2, ub = ub_fac, nev = 3, seed = 0);
            push!(spec_cut_list, spec_cut)
            push!(hmetis_cut_list, hmetis_cut)
        end

        local df = DataFrame(:Testcase => test_case_list, :hMetis_Cut => hmetis_cut_list, :Refined_Cut => spec_cut_list)
        CSV.write(csv_name, df)
    end
end=#


files = readdir("/home/bodhi91/DEV/sandbox/Partitioner/SmartPartAndExRefine/ISPD_benchmark/hmetis_format/", join=true)
file_names = readdir("/home/bodhi91/DEV/sandbox/Partitioner/SmartPartAndExRefine/ISPD_benchmark/hmetis_format/")
nparts_list = [2]
ub_list = [2, 5, 10]
nev_list = [3]

t = @elapsed for ub_fac in ub_list
    for n_eig in nev_list
        csv_name = "ISPD_" * string(ub_fac) * "_" * string(n_eig) * "_results.csv" 
        test_case_list = String[]
        spec_cut_list = Int[]
        hmetis_cut_list = Int[]

        for i in 1:length(files)
            partition_vector = Int[]
            f = files[i]
            push!(test_case_list, file_names[i])
            @info "Running $f with ub_fac = $ub_fac, eigen_vecs = $n_eig"
            ~, spec_cut, hmetis_cut  = Main.SpectralRefinement.SpectralHmetisRefinement(hg = f, nparts = 2, ub = ub_fac, nev = 3, seed = 0);
            push!(spec_cut_list, spec_cut)
            push!(hmetis_cut_list, hmetis_cut)
        end

        local df = DataFrame(:Testcase => test_case_list, :hMetis_Cut => hmetis_cut_list, :Refined_Cut => spec_cut_list)
        CSV.write(csv_name, df)
    end
end


#=t = @elapsed for ub_fac in ub_list
    csv_name = string(ub_fac) * "_" * "hmetis_results.csv" 
    test_case_list = String[]
    hmetis_cut_list = Int[]

    for i in 1:length(files)
        partition_vector = Int[]
        f = files[i]
        push!(test_case_list, file_names[i])
        @info "Running $f with ub_fac = $ub_fac"
        (hedges, eptr, vertex_weights, hyperedge_weights, num_vertices, num_hyperedges, ~) = Main.MultilevelPartitioner.ReadHypergraphFile(f)
        hypergraph = Main.MultilevelPartitioner.Hypergraph(num_vertices, num_hyperedges, hedges, eptr, vertex_weights, hyperedge_weights)
        Main.MultilevelPartitioner.hmetis(f, 2, ub_fac, 10, 1, 1, 0, 1, 0)
        fpart_name = f * ".part.2"
        fp = open(fpart_name)

        for ln in eachline(fp)
            i += 1
            push!(partition_vector, parse(Int, ln))
        end

        close(fp)

        cut_size, ~ = Main.MultilevelPartitioner.FindCutSize(partition_vector, hypergraph)
        @info "Cut size: $cut_size"

        push!(hmetis_cut_list, cut_size)

        cmd = "rm " * fpart_name
        run(`sh -c $cmd`)
    end

    local df = DataFrame(:Testcase => test_case_list, :Hmetis_Cut => hmetis_cut_list)
    CSV.write(csv_name, df)
end=#

#=
files = readdir("/home/bodhi91/DEV/sandbox/Partitioner/Benchmarks/XLNX_Benchmarks_hMetis_format/", join=true)
file_names = readdir("/home/bodhi91/DEV/sandbox/Partitioner/Benchmarks/XLNX_Benchmarks_hMetis_format/")

t = @elapsed for ub_fac in ub_list
    for n_eig in nev_list
        csv_name = "XLNX_" * string(ub_fac) * "_" * string(n_eig) * "_results.csv" 
        test_case_list = String[]
        spec_cut_list = Int[]

        for i in 1:length(files)
            partition_vector = Int[]
            f = files[i]
            push!(test_case_list, file_names[i])
            @info "Running $f with ub_fac = $ub_fac, eigen_vecs = $n_eig"
            part_vec, cut_size = Main.MultilevelPartitioner.SpectralPart(hg = f, nparts = 2, ub = ub_fac, nev = n_eig, seed = 0)
            push!(spec_cut_list, cut_size)
        end

        local df = DataFrame(:Testcase => test_case_list, :Spectral_Cut => spec_cut_list)
        CSV.write(csv_name, df)
    end
end

@info "Total Time for Massive Exp: $t"


# hmetis default run set up for 40 starts

ub_factor = 2

for f in files
    for i in 1:40
        hmetis(f, 2, ub_factor, 10, 1, 1, 0, 1, 0)
        
        cmd = "mv " * f * ".part.2" * string(i) * "_" * f * ".part.2"
        run(`sh -c $cmd`, wait=true)
    end

end
=#