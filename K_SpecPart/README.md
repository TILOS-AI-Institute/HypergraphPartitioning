# K-SpecPart #

K-SPecPart is an extension (and improvement) over SpecPart that can compute multi-way (K >= 2) partitions. 

# Requirements #

In order to run K-SpecPart, the following requirements are mandatory: 
* [Julia](https://julialang.org/) version 1.6.0 or newer. 
* The following Julia libraries:
  * [Shuffle](https://docs.juliahub.com/Shuffle/X0eqg/0.1.1/)
  * [LightGraphs](https://github.com/sbromberger/LightGraphs.jl)
  * [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl)
  * [LDLFactorizations](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl)
  * [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/)
  * [Laplacians](https://github.com/danspielman/Laplacians.jl)
  * [LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl)
  * [MultivariateStats](https://github.com/JuliaStats/MultivariateStats.jl)
* [OR-Tools](https://developers.google.com/optimization)
* If you have acces to [CPLEX ILP solver](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290), you can use CPLEX in place of OR-Tools. This change can be made by modifications to ```CMakeLists.txt``` in ```ilp_partitioner``` directory. 
* [METIS](https://github.com/KarypisLab/METIS)
* [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)
* [FM](https://github.com/ABKGroup/TritonPart_OpenROAD.git)

Julia libraries can be installed by doing the following in Julia REPL:

```
julia> using Pkg
julia> Pkg.add("Package name")
```

# K-SpecPart Parameters #

K-SpecPart accepts the following parameters:

| Parameter   | Description |
| ----------- | ----------- |
| hypergraph_fixed_file      | The fixed file consisting of fixed vertices       |
| hint_file   | The hint partition file        |
| imb  | The imbalance factor         |
| num_parts   | The number of partitions        |
| eigvecs   | The number of eigenvectors (Default 2)        |
| refine_iters   | The number of K-SpecPart iterations (Default 2)        |
| solver_iters   | The number of LOCPCG solver iterations (Default 40)        |
| best_solns   | The number of partition solutions picked for overlay-clustering (Default 5)        |
| ncycles   | The number of random cycles for constructing spectral sparsifier graph (Default 2)        |
| seed   | The random seed for K-SpecPart (Default 0)       |

# K-SpecPart Runtime and cutsize comparison for K = 2 #
| Benchmark     | Vertices | Hyperedges | Avg hMETIS time (s) | K-SpecPart time (s) | Avg hMETIS cutsize | K-SpecPart cutsize |
|---------------|----------|------------|-----------------|-------------------|----------------|------------------| 
| sparcT1_core  | 91976    | 92827      | 9.5            | 32.52             | 982.0           | 979             |
| neuron        | 92290    | 125305     | 6.3            | 33.14              | 245.0            | 244              |
| stereovision  | 94050    | 127085     | 8.3            | 22.19             | 171.0            | 169              | 
| des90         | 111221   | 139557     | 10.2             | 18.22             | 377.5            | 380              | 
| SLAM_spheric  | 113115   | 142408     | 12.7           | 45.18             | 1061.0           | 1061             | 
| cholesky_mc   | 113250   | 144948     | 8.6            | 30.36             | 282.5            | 282              | 
| segmentation  | 138295   | 179051     | 13.2            | 47.12             | 120.2            | 120              | 
| bitonic_mesh  | 192064   | 235328     | 17.6           | 64.24             | 585.2            | 584              | 
| dart          | 202354   | 223301     | 14.4           | 51.22             | 841.3            | 793              | 
| openCV        | 217453   | 284108     | 11.7           | 42.17             | 435.3            | 434              | 
| stap_qrd      | 240240   | 290123     | 15.9           | 56.31             | 378.0            | 373              | 
| minres        | 261359   | 320540     | 18.2           | 82.92             | 207.0            | 207              | 
| cholesky_bdti | 266422   | 342688     | 19.8           | 91.33             | 1156.0           | 1136             | 
| denoise       | 275638   | 356848     | 26.3           | 95.75             | 477.0            | 418              | 
| sparcT2_core  | 300109   | 302663     | 29.7           | 112.31             | 1216.3           | 1193             |
| gsm_switch    | 493260   | 507821     | 27.2           | 110.25             | 4210.3           | 1835             | 
| mes_noc       | 547544   | 577664     | 42.3           | 121.92             | 634.0            | 634              | 
| LU230         | 574372   | 669477     | 46.6           | 162.32             | 3328.0           | 3263             | 
| LU_Network    | 635456   | 726999     | 64.1           | 182.23             | 524.0            | 524              | 
| sparcT1_chip2 | 820886   | 821274     | 76.8           | 216.82            | 894.7           | 876              | 
| directrf      | 931275   | 1374742    | 85.6           | 223.47            | 630.6            | 574              | 
| bitcoin_miner | 1089284  | 1448151    | 89.7           | 369.19            | 1299.3           | 1297             | 


# K-SpecPart Runtime and cutsize comparison for K = 3 #
| Benchmark     | Vertices | Hyperedges | Avg hMETIS time (s) | K-SpecPart time (s) | Avg hMETIS cutsize | K-SpecPart cutsize |
|---------------|----------|------------|-----------------|-------------------|----------------|------------------| 
| sparcT1_core  | 91976    | 92827      | 12.2            | 60.39             | 2530.0           | 1878             |
| neuron        | 92290    | 125305     | 16.7            | 70.55              | 366.3            | 396              |
| stereovision  | 94050    | 127085     | 14.1            | 64.06             | 343.0            | 344              | 
| des90         | 111221   | 139557     | 14.8             | 79.95             | 524.0            | 535              | 
| SLAM_spheric  | 113115   | 142408     | 19.2           | 72.23             | 2724.5           | 2720             | 
| cholesky_mc   | 113250   | 144948     | 21.3            | 82.79             | 867.7            | 859              | 
| segmentation  | 138295   | 179051     | 22.4            | 86.86             | 447.0            | 438              | 
| bitonic_mesh  | 192064   | 235328     | 23.7           | 107.24             | 895.0            | 895              | 
| dart          | 202354   | 223301     | 29.1           | 112.64             | 1226.9            | 1243              | 
| openCV        | 217453   | 284108     | 28.8           | 137.72             | 525.0            | 525              | 
| stap_qrd      | 240240   | 290123     | 39.3           | 153.39             | 505.7            | 501              | 
| minres        | 261359   | 320540     | 42.4           | 167.89             | 309.0            | 309              | 
| cholesky_bdti | 266422   | 342688     | 41.8           | 171.23             | 1769.2           | 1700             | 
| denoise       | 275638   | 356848     | 37.9           | 189.67             | 888.3            | 918              | 
| sparcT2_core  | 300109   | 302663     | 46.7           | 176.26             | 2788.9           | 2521             |
| gsm_switch    | 493260   | 507821     | 52.1           | 172.92             | 4328.2           | 3702             | 
| mes_noc       | 547544   | 577664     | 49.8           | 181.42             | 1164.2            | 1167              | 
| LU230         | 574372   | 669477     | 48.6           | 192.18             | 4572.6           | 4570             | 
| LU_Network    | 635456   | 726999     | 52.3           | 210.04             | 882.0            | 882              | 
| sparcT1_chip2 | 820886   | 821274     | 55.8           | 212.23            | 1418.7           | 1365              | 
| directrf      | 931275   | 1374742    | 59.2           | 222.18            | 759.2            | 761              | 
| bitcoin_miner | 1089284  | 1448151    | 60.3           | 265.86            | 1964.3           | 1917             | 


# K-SpecPart Runtime and cutsize comparison for K = 4 #
| Benchmark     | Vertices | Hyperedges | Avg hMETIS time (s) | K-SpecPart time (s) | Avg hMETIS cutsize | K-SpecPart cutsize |
|---------------|----------|------------|-----------------|-------------------|----------------|------------------| 
| sparcT1_core  | 91976    | 92827      | 29.1            | 74.41             | 2543.3           | 2449             |
| neuron        | 92290    | 125305     | 23.7            | 84.64              | 435.5            | 429              |
| stereovision  | 94050    | 127085     | 21.5            | 79.42             | 471.2            | 476              | 
| des90         | 111221   | 139557     | 20.8             | 89.25             | 725.4            | 747              | 
| SLAM_spheric  | 113115   | 142408     | 22.9           | 97.31             | 3336.7           | 3274             | 
| cholesky_mc   | 113250   | 144948     | 21.2            | 101.24             | 981.1            | 984              | 
| segmentation  | 138295   | 179051     | 33.2            | 116.52             | 495.7            | 493              | 
| bitonic_mesh  | 192064   | 235328     | 37.8           | 130.87             | 1304.7            | 1311              | 
| dart          | 202354   | 223301     | 47.3           | 145.51             | 1456.2            | 1455              | 
| openCV        | 217453   | 284108     | 52.9           | 192.14             | 526.2            | 521              | 
| stap_qrd      | 240240   | 290123     | 55.2           | 207.81             | 707.7            | 641              | 
| minres        | 261359   | 320540     | 56.9           | 184.54             | 407.0            | 405              | 
| cholesky_bdti | 266422   | 342688     | 59.1           | 174.65             | 1870.3           | 1869             | 
| denoise       | 275638   | 356848     | 66.3           | 182.19             | 1199.7            | 1186              | 
| sparcT2_core  | 300109   | 302663     | 69.1           | 198.42             | 3592.9           | 3549             |
| gsm_switch    | 493260   | 507821     | 67.2           | 202.65             | 5035.2           | 4049             | 
| mes_noc       | 547544   | 577664     | 62.7           | 201.18             | 1310.3            | 1357              | 
| LU230         | 574372   | 669477     | 66.9           | 204.73             | 6320.1           | 6311             | 
| LU_Network    | 635456   | 726999     | 71.2           | 265.67             | 1498.4            | 1417              | 
| sparcT1_chip2 | 820886   | 821274     | 72.8           | 395.13            | 1764.2           | 1757              | 
| directrf      | 931275   | 1374742    | 85.6           | 301.45            | 1216.4            | 793              | 
| bitcoin_miner | 1089284  | 1448151    | 89.7           | 312.23            | 1962.7           | 1859             | 

