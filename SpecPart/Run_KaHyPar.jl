files = readdir("/home/bodhi91/sandbox/Partitioners/Benchmarks/Xilinx", join=true)
eps = [0.02, 0.04, 0.1, 0.2]
config = "cut_rKaHyPar_sea20.ini"

for hypergraph in files
    for ep in eps
        @info "Running $hypergraph with eps $ep"
        cmd = "./KaHyPar -h " * hypergraph * " -k 2 -e " * string(ep) * " -o cut -m recursive -w true -p " * config
        run(`sh -c $cmd`, wait=true)
    end
end
