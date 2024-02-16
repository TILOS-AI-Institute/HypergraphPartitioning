# SpecPart User Guide #

SpecPart is the first supervised spectral hypergraph partitioning framework which enhances a given hypergraph partition solution.
This directory contains the Julia implementation of SpecPart. SpecPart can be accessed [here](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/SpecPart/SpecPart_final_submission.pdf). 

# Requirements #

In order to run SpecPart, the following requirements are mandatory: 
* [Julia](https://julialang.org/) version 1.6.0 or newer. 
* The following Julia libraries:
  * [Shuffle](https://docs.juliahub.com/Shuffle/X0eqg/0.1.1/)
  * [LightGraphs](https://github.com/sbromberger/LightGraphs.jl)
  * [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl)
  * [LDLFactorizations](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl)
  * [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/)
  * [Laplacians](https://github.com/danspielman/Laplacians.jl)
  * [LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl)
* [CPLEX ILP solver](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290) 
* [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)

Julia libraries can be installed by doing the following in Julia REPL:

```
julia> using Pkg
julia> Pkg.add("Package name")
```

# SpecPart Parameters #

SpecPart accepts the following parameters:

| Parameter   | Description |
| ----------- | ----------- |
| hg      | The hypergraph file.       |
| pfile   | The partition file.        |
| Nparts   | The number of partitions. (Current implementation only allows bipartitions)         |
| cycles   | The number of random cycles used for graph construction. (Default 2)        |
| hyperedges_threshold   | The threshold number of hyperedges used for running ILP. (Default 300)        |
| ub   | The imbalance factor.        |
| nev   | The number of eigenvectors to be computed. (Default 2)        |
| refine_iters   | The number of ISSHP iterations. (Default 2)        |
| best_solns   | The number of partition solutions picked for overlay-clustering. (Default 5)        |
| seed   | The random seed for SpecPart. (Default 0)       |
| solver_iters   | The number of iterations of LOBPCG. (Default 80)        |

# SpecPart Example #
We show how we run SpecPart on ISPD98 testcase "ibm02.hgr" with an input partition. As seen from this example SpecPart improves the initial partition of 339 to 336.

```
julia> include("SpecPart/SpectralRefinement.jl")
julia> using Main.SpecPart
julia> SpecPart.SpectralHmetisRefinement(hg = "benchmark/ISPD_benchmark/ibm02.hgr", pfile = "solutions/ISPD_benchmark_solutions/hMetis/UBfactor_2/ibm02.hgr.k.2.UBfactor.2.seed.0", 
Nparts = 2, cycles = 2, hyperedges_threshold = 300, ub = 2, nev = 2, refine_iters = 2, best_solns = 5)

[ Info: ================================================================================
[ Info: STARTING SPECPART: SUPERVISED SPECTRAL PARTITIONING ENGINE
[ Info: ================================================================================
[ Info: TOTAL VERTICES: 19601
[ Info: TOTAL HYPEREDGES: 19584
[ Info: POST PROCESSING :: TOTAL VERTICES: 19601
[ Info: POST PROCESSING :: TOTAL HYPEREDGES: 19584
[ Info: SIZE OF LARGEST HYPEREDGE: 134
[ Info: MAX CAPACITY CONSTRAINT: 10193
[ Info: MIN CAPACITY CONSTRAINT: 9408
[ Info: ================================================================================
[ Info: [HINT] CUT RECORDED IS 339
┌ Info: [EIGEN VECTOR DETAILS] :: Results of LOBPCG Algorithm
│  * Algorithm: LOBPCG - CholQR
│  * λ: [6.09867746662552e-8,8.185494183481112e-8]
│  * Residual norm(s): [1.3833209684774178e-17,1.3317377476442714e-16]
│  * Convergence
│    * Iterations: 81
│    * Converged: false
└    * Iterations limit: 80
[ Info: [EIGEN VECTOR DETAILS] :: SOLVER TIME :: 1.712944723 seconds
[ Info: [SPECTRAL ITERATION 1] MIN CUT FOUND :: 434 :: TREE GEN-SOLVE TIME :: 1.092160598 seconds
┌ Info: [EIGEN VECTOR DETAILS] :: Results of LOBPCG Algorithm
│  * Algorithm: LOBPCG - CholQR
│  * λ: [6.063794065520354e-8,8.233442427641748e-8]
│  * Residual norm(s): [9.10638903776076e-16,1.0476424071977864e-14]
│  * Convergence
│    * Iterations: 81
│    * Converged: false
└    * Iterations limit: 80
[ Info: [EIGEN VECTOR DETAILS] :: SOLVER TIME :: 1.905182079 seconds
[ Info: [SPECTRAL ITERATION 2] MIN CUT FOUND :: 431 :: TREE GEN-SOLVE TIME :: 0.68814692 seconds
[ Info: ================================================================================
[ Info: SIZE OF CLUSTERED HYPERGRAPH IS: 952 VERTICES AND 875 HYPEREDGES
[ Info: ================================================================================
[ Info: PARTITIONING CLUSTERED HYPERGRAPH .....
[ Info: RUNNING CPLEX AS GOLDEN PARTITIONER
[ Info: [POST TOOL CUT] CUT RECORDED IS 336
[ Info: ================================================================================
[ Info: [SUPERVISED SPECTRAL] CUTSIZE OBTAINED: 336
[ Info: [SUPERVISED SPECTRAL] AREA SPLIT OBTAINED: 10191 [0.52%] :: 9410 [0.48%]
[ Info: ================================================================================
[ Info: [RUNTIME] IO PROCESSING :: 0.11667171 SECONDS
[ Info: [RUNTIME] SPECTRAL :: 5.562175724 SECONDS
[ Info: [RUNTIME] CLUSTERING :: 0.048536965 SECONDS
[ Info: [RUNTIME] PARTITION CLUSTERED HYPERGRAPH :: 1.311114935 SECONDS
[ Info: ================================================================================
[ Info: [RUNTIME] TOTAL EXECUTION TIME :: 7.038499334000001 SECONDS
[ Info: ================================================================================

```

# Experimental Results #

Runtime comparison of SpecPart vs hMETIS on ISPD98 benchmark suite with unit vertex weights and actual vertex weights for an imbalance factor of 2: 

|           |          |            | Unit weight     | Unit weight       | Actual weight   | Actual weight     |
|-----------|----------|------------|-----------------|-------------------|-----------------|-------------------|
| Benchmark | Vertices | Hyperedges | hMETIS time (s) | SpecPart time (s) | hMETIS time (s) | SpecPart time (s) |
| IBM01     | 12752    | 14111      | 1.31            | 2.6               | 1.35            | 1.46              |
| IBM02     | 19601    | 19584      | 2.54            | 5.65              | 2.49            | 5.63              |
| IBM03     | 23136    | 27401      | 2.56            | 5.92              | 2.61            | 5.73              |
| IBM04     | 27507    | 31970      | 2.53            | 6.32              | 3.02            | 7.26              |
| IBM05     | 29347    | 28446      | 3.77            | 5.98              | 3.72            | 8.12              |
| IBM06     | 32498    | 34826      | 3.15            | 6.17              | 4.08            | 7.82              |
| IBM07     | 45926    | 48117      | 4.93            | 12.08             | 6.14            | 10.46             |
| IBM08     | 51309    | 50513      | 6.67            | 13.48             | 7.86            | 12.58             |
| IBM09     | 53395    | 60902      | 4.5             | 8.83              | 5.68            | 13.42             |
| IBM10     | 69429    | 75196      | 8.3             | 14.24             | 9.23            | 16.53             |
| IBM11     | 70558    | 81454      | 7.15            | 14.72             | 7.94            | 15.02             |
| IBM12     | 71076    | 77240      | 9.18            | 18.68             | 10.28           | 19.74             |
| IBM13     | 84199    | 99666      | 7.74            | 14.42             | 8.79            | 16.44             |
| IBM14     | 147605   | 152772     | 17.62           | 30.32             | 19.59           | 37.73             |
| IBM15     | 161570   | 186608     | 21.05           | 44.68             | 22.18           | 44.31             |
| IBM16     | 183484   | 190048     | 22.27           | 48.19             | 26.51           | 47.83             |
| IBM17     | 185495   | 189581     | 30.49           | 52.01             | 35.16           | 49.76             |
| IBM18     | 210613   | 201920     | 27.54           | 52.84             | 29.58           | 54.27             |


Runtime and cutsize comparison of SpecPart vs hMETIS on Titan23 benchmark suite for an imbalance factor of 2: 

| Benchmark     | Vertices | Hyperedges | hMETIS time (s) | SpecPart time (s) | hMETIS cutsize | SpecPart cutsize | % Improvement |
|---------------|----------|------------|-----------------|-------------------|----------------|------------------| --------------|
| sparcT1_core  | 91976    | 92827      | 7.05            | 19.61             | 1073           | 1012             | 5.68          |
| neuron        | 92290    | 125305     | 4.75            | 9.26              | 260            | 252              | 3.08          |
| stereovision  | 94050    | 127085     | 5.23            | 15.52             | 213            | 180              | 15.49         |
| des90         | 111221   | 139557     | 9.2             | 18.22             | 403            | 402              | 0.25          |
| SLAM_spheric  | 113115   | 142408     | 10.23           | 25.98             | 1061           | 1061             | 0             |
| cholesky_mc   | 113250   | 144948     | 6.11            | 10.64             | 301            | 285              | 5.31          |
| segmentation  | 138295   | 179051     | 9.77            | 16.23             | 141            | 126              | 10.63         |
| bitonic_mesh  | 192064   | 235328     | 15.05           | 31.46             | 667            | 585              | 12.29         |
| dart          | 202354   | 223301     | 13.38           | 26.82             | 849            | 807              | 4.95          |
| openCV        | 217453   | 284108     | 10.68           | 22.39             | 535            | 510              | 4.67          |
| stap_qrd      | 240240   | 290123     | 14.52           | 31.98             | 399            | 399              | 0             |
| minres        | 261359   | 320540     | 17.01           | 30.28             | 215            | 215              | 0             |
| cholesky_bdti | 266422   | 342688     | 18.17           | 36.42             | 1161           | 1156             | 0.43          |
| denoise       | 275638   | 356848     | 23.24           | 44.03             | 814            | 416              | 58.89         |
| sparcT2_core  | 300109   | 302663     | 25.59           | 51.28             | 1282           | 1244             | 2.96          |
| gsm_switch    | 493260   | 507821     | 21.67           | 43.87             | 5883           | 1852             | 68.52         |
| mes_noc       | 547544   | 577664     | 36.82           | 64.83             | 674            | 641              | 4.89          |
| LU230         | 574372   | 669477     | 37.46           | 68.12             | 3328           | 3273             | 1.61          |
| LU_Network    | 635456   | 726999     | 57.06           | 71.24             | 549            | 525              | 4.37          |
| sparcT1_chip2 | 820886   | 821274     | 57.06           | 123.12            | 1198           | 899              | 24.92         |
| directrf      | 931275   | 1374742    | 71.34           | 143.19            | 588            | 574              | 2.38          |
| bitcoin_miner | 1089284  | 1448151    | 83.26           | 165.34            | 1576           | 1514             | 3.93          |
 

# Reference #

To reference this work, please cite: 

I. Bustany, A. B. Kahng, Y. Koutis, B. Pramanik and Z. Wang, "SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement", Proc. ACM/IEEE International Conference on Computer-Aided Design, 2022. (https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/SpecPart/SpecPart_final_submission.pdf)
