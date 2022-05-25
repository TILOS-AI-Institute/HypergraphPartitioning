using DataFrames, CSV

include("SpectralRefinement.jl")

files = readdir("/home/bodhi91/DEV/sandbox/Partitioner/Benchmarks_Spectral/titan23_benchmark/CUHK_benchmark/", join=true)
file_names = readdir("/home/bodhi91/DEV/sandbox/Partitioner/Benchmarks_Spectral/titan23_benchmark/CUHK_benchmark/")
nparts_list = [2]
ub_list = [1, 2, 5, 10]
nev_list = [3]

#files = readdir("../SmartPartAndExRefine/ISPD_benchmark/hmetis_format/", join=true)

t = @elapsed for ub_fac in ub_list
    for n_eig in nev_list
        csv_name = "ex_refined_TITAN_" * string(ub_fac) * "_" * string(n_eig) * "_results.csv" 
        test_case_list = String[]
        spec_cut_list = Int[]
        hmetis_cut_list = Int[]

        for i in 1:length(files)
            f = files[i]
            pname = "/home/bodhi91/DEV/sandbox/Partitioner/IterativeSpectralRefinement/PartitionCertificates/UB_factor_" * string(ub_fac) * "/" * file_names[i] *".part.2"
            push!(test_case_list, file_names[i])
            @info "Running $f with ub_fac = $ub_fac, eigen_vecs = $n_eig"
            ~, spec_cut, hmetis_cut = Main.SpectralRefinement.SpectralHmetisRefinement(hg = f, pfile = pname, nparts = 2, ub = ub_fac, nev = 3, seed = 0);
            @info "Spectral cut :: $spec_cut       Hmetis cut :: $hmetis_cut"
            push!(spec_cut_list, spec_cut)
            push!(hmetis_cut_list, hmetis_cut)
        end

        local df = DataFrame(:Testcase => test_case_list, :hMetis_Cut => hmetis_cut_list, :Refined_Cut => spec_cut_list)
        CSV.write(csv_name, df)
    end
end